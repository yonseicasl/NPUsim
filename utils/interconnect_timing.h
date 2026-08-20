#ifndef __INTERCONNECT_TIMING_H__
#define __INTERCONNECT_TIMING_H__

#include <cstddef>
#include <vector>
#include "def.h"

// Timing-only contracts shared by interconnect hierarchy components.
// STORE_AND_FORWARD is modeled as a serialized unicast fabric until a routed
// mesh model is introduced; unsupported topologies must fail at initialization.
bool is_supported_spatial_noc(noc_type_t topology);
bool has_valid_active_shape(unsigned physical_height, unsigned physical_width,
                            unsigned active_height, unsigned active_width);
bool is_valid_memory_line_bits(size_t width_bits);

// L2 (CE5/P1-C): PHYSICAL packet decomposition of one source->link->destination transfer.
//
// A width-conversion boundary is a barrier only between the transactions that share a
// packet: with an 8b source line feeding a 256b link, 32 source reads must land before
// THAT link transaction fires, but the next 32 reads may proceed while it is in flight.
// The packet is therefore the WIDEST granularity on the path (source line, link, or
// destination line): one transaction of the widest stage spans a whole number of every
// narrower stage's transactions, which is exactly where packing/unpacking happens, and
// the widest stage contributes one transaction per packet.
//
// Packets come from the payload and the widths -- NOT from an even split of the aggregate
// per-stage counts. An even split invents boundaries no line width produces: 200 bytes
// over an 8b line into 256b link packets is 32,32,32,32,32,32,8, never 29,29,29,29,28,28,28.
// All packets but the last are full; the last carries the remainder.
struct transfer_packet_groups_t {
    size_t packets;                   // full packets plus the tail packet (0 when nothing moves)
    size_t source_per_packet;         // counts in a full packet
    size_t link_per_packet;
    size_t destination_per_packet;
    size_t source_tail;               // counts in the last packet (== *_per_packet when exact)
    size_t link_tail;
    size_t destination_tail;

    size_t source_transactions() const;
    size_t link_transactions() const;
    size_t destination_transactions() const;
    // A physically absent stage is dropped from the pipeline, not charged as zero-cost
    // work: a pass-through PE array has no temporal-buffer write, and a GLB-bypassed
    // stream has no GLB read.
    transfer_packet_groups_t without_source() const;
    transfer_packet_groups_t without_destination() const;
    // Replace the destination stage's TOTAL transaction count, keeping the packet
    // boundaries of the shared source/link payload. For a destination that moves its own
    // SHARE of the payload rather than all of it (PE-array distribution: each PE writes
    // only its slice into its own local buffer), so the exact per-packet destination
    // grouping is not derivable from the shared payload; the total is spread over the
    // packets instead, larger packets first.
    transfer_packet_groups_t with_destination_total(size_t transactions) const;
};

// `payload_bits` is the transfer's storage size; the three widths are the source line,
// the link, and the destination line. The aggregate counts are the ones the caller charges
// cycles and energy from: the decomposition is reconciled against them so the makespan can
// never describe more or less work than was billed (any residue -- a separately counted
// metadata stream, or widths that are not integer multiples of each other -- lands in the
// tail packet).
transfer_packet_groups_t transfer_packet_groups(size_t payload_bits,
                                                size_t source_line_bits, size_t link_bits,
                                                size_t destination_line_bits,
                                                size_t source_transactions,
                                                size_t link_transactions,
                                                size_t destination_transactions);

// Makespan of the decomposed transfer, as a resource-availability recurrence over the
// packets. Each stage is one serial resource; within a packet a stage may start as soon as
// (a) its own resource is free and (b) its dependency is met -- the WHOLE upstream packet
// for a gather boundary (fewer downstream transactions than upstream), the FIRST upstream
// transaction for a fan-out boundary -- and its last transaction can never precede the
// last upstream transaction that feeds it. Consecutive packets therefore overlap
// stage-wise, which a "first packet latency + sum of later packet periods" formula gets
// wrong as soon as the packets differ in size (a cheap tail packet must not be charged a
// full packet's period, and it can still be gated by a downstream stage that is behind).
double pipelined_transfer_cycles(const transfer_packet_groups_t &groups,
                                 double source_stage, double link_stage,
                                 double destination_stage);
// Equal-width endpoints: `transactions` packets of one transaction per stage.
double pipelined_transfer_cycles(unsigned transactions, double source_stage,
                                 double link_stage, double destination_stage);

// P3: streaming-pipeline makespan of a set of stage-busy TOTALS that overlap over
// the SAME `tiles` temporal repetitions (e.g. GLB/PE-array/PE all continuously
// processing the same GLB-repetition count). PRECONDITION: every value in
// stage_totals must scale at the SAME repetition rate as `tiles` -- do not mix a
// stage whose total scales per-datatype (e.g. DRAM/multi-chip access cycles, which
// scale with a per-datatype repetition factor, not the uniform GLB-repetition count)
// into one call with stages that scale uniformly by `tiles`, or the derived per-tile
// costs are meaningless. With <=1 tile there is nothing to pipeline against (no tile
// i+1 for an upstream stage to work on while a downstream stage handles tile i), so
// this degenerates to the existing flat max() -- NOT merely a conservative choice:
// the general formula collapses to a flat SUM at tiles<=1, which would silently
// contradict an "these stages overlap" boundary.
//
// L5 NOTE: superseded for the layer critical path by pipeline_timeline_cycles() below,
// which is a real per-tile event recurrence instead of an average-per-tile closed form.
// Kept for the fill+bottleneck contract it documents and for its unit tests.
double temporal_pipeline_run_cycles(unsigned tiles, const std::vector<double> &stage_totals);

// L5: per-tile pipeline timeline over a linear chain of stages with BOUNDED buffers.
//
// `stage_tile_costs[s][k]` is stage s's cost for tile k -- a per-tile cost VECTOR, not a
// single average, so a stage that serves only some tiles (an off-chip stage refetching at a
// lower rate than the on-chip stages consume) carries 0 on the tiles it does not serve, and
// a cheaper tail tile is expressible where the mapping produces one. Every stage must have
// the same vector length: the tile index is a shared, normalized clock.
//
// `boundary_depths[b]` is how many tiles may be in flight across the boundary from stage b
// to stage b+1 -- i.e. the staging buffer's depth in tiles:
//   1 : single buffer. The producer cannot run ahead: tile k+1 may not start until the
//       consumer has taken tile k. This is the "no overlap at this boundary" case.
//   2 : double buffer. The producer fills one half while the consumer drains the other --
//       the classic overlap.
//   N : an N-deep queue; the producer may run N-1 tiles ahead.
// A depth of 0 is treated as 1 (there is always at least one tile in flight).
//
// The recurrence, per tile k and stage s:
//   start(s,k)  = max( finish(s,k-1),                 // the stage is one serial resource
//                      finish(s-1,k),                 // its input tile must have arrived
//                      finish(s+1,k-depth[s]) )       // and a downstream slot must be free
//   finish(s,k) = start(s,k) + cost(s,k)
// The third term is the back-pressure that an average-per-tile closed form cannot express:
// a fast producer in front of a slow consumer stalls once the queue fills, rather than
// running to completion in parallel.
//
// `stall` (optional, one entry per stage) receives the time each stage spent blocked
// specifically by that back-pressure term, so a report can say WHERE the pipeline stalled
// rather than only how long it took.
double pipeline_timeline_cycles(const std::vector<std::vector<double>> &stage_tile_costs,
                                const std::vector<unsigned> &boundary_depths,
                                std::vector<double> *stall = 0);

// CE4/P1-D: combine a per-datatype cost vector for ONE entity (buffer instance / port).
// `serialized_types` selects the rule: a shared port serializes the three datatype
// streams (sum), separate per-datatype partitions serve them concurrently (max).
double combine_datatype_cycles(const std::vector<double> &per_type, bool serialized_types);

// CE4/P1-D: reduce a per-entity, per-datatype cost matrix into one stage axis in the ONLY
// correct order -- combine each entity's OWN datatypes first, then take the max ACROSS
// entities: max_entity(combine_type(per_entity_type[e][t])). Collapsing entities first
// (per_type[t] = max_entity(...), then combining types) yields combine_type(max_entity(...)),
// which invents elapsed time that never existed: with a shared GLB where chip 0 spends 100
// cycles filling inputs and chip 1 spends 100 filling weights, the wrong order reports 200
// while the two chips actually run in parallel and finish in 100.
double entity_combined_cycles(const std::vector<std::vector<double>> &per_entity_type,
                              bool serialized_types);

// L1: same reduction for a port with TWO cost sides that must be added per datatype --
// e.g. a global buffer's read side (GLB -> PE array) and its fill side (multi-chip -> GLB),
// which are tracked apart only because they scale by different repetition factors:
//     max_entity(combine_type(a[e][t] + b[e][t]))
// Adding the already-combined sides instead (combine_type(a) + combine_type(b)) is WRONG on
// a SEPARATE buffer: max_type(a) and max_type(b) can be peaks of DIFFERENT datatype
// partitions, and those partitions run in parallel, so summing them charges elapsed time
// that never happened. (On a shared buffer both combinations are sums and the two forms
// agree, which is why this stayed hidden.) A shorter/absent matrix contributes zeros.
double entity_combined_cycles(const std::vector<std::vector<double>> &a,
                              const std::vector<std::vector<double>> &b,
                              bool serialized_types);

// L8: analytical open-page row-activation cost of ONE sequential dense stream of
// `stream_bytes` bytes over `row_buffer_bytes`-byte rows spread across `num_banks` banks.
//
// `activations` is how many rows the stream opens, and energy scales with it -- every
// activation costs energy whether or not it overlapped another bank's. `busiest_bank` is how
// many of those activations the busiest bank has to serialize, and latency scales with THAT.
//
// The bank term is an IDEAL even spread, not an address map: this model has no per-request
// address, so it cannot say which bank a real access would conflict on. It is therefore a
// lower bound on activation latency, and the report states that scope explicitly (see
// dram_t::describe_timing_limits()). A zero row size disables the model.
struct dram_row_activation_cost_t {
    size_t activations;
    size_t busiest_bank;
};

dram_row_activation_cost_t dram_row_activations(size_t stream_bytes, size_t row_buffer_bytes,
                                                unsigned num_banks);

// Unit static energy is expressed in pJ/cycle.
double static_energy_for_cycles(double unit_energy, double elapsed_cycles);

// All widths below are expressed in bits.  The payload/metadata counters are
// reported independently even for interleaved metadata layouts; `link` is the
// physical serialized transaction count used for timing and energy.
struct datatype_transfer_timing_t {
    size_t source_accesses;
    size_t destination_accesses;
    size_t link_transactions;
    size_t pipeline_transactions;
    size_t payload_link_transactions;
    size_t metadata_link_transactions;
    // L2: the physical packet decomposition of this transfer, derived from the same
    // payload and widths as the counts above. Feed this to pipelined_transfer_cycles()
    // instead of the raw counts, which cannot express where the packet boundaries are.
    transfer_packet_groups_t groups;
};

datatype_transfer_timing_t datatype_transfer_timing(data_type_t type, size_t elements,
                                                     size_t source_line_bits,
                                                     size_t destination_line_bits,
                                                     size_t link_bits);
// Pipelined-hop contract (SP1): a stream of T transactions to the farthest active
// destination completes in (T + latency_fill_hops)*noc_cycle -- hops pipeline, so
// the route depth is a one-time fill, not a per-transaction multiplier.
//
// E2: LATENCY and ENERGY need different quantities, and the two DIRECTIONS of traffic across a
// spatial array need different ones again. Energy is per physical link crossed, so it must
// count the links a delivery actually uses:
//
//   MULTICAST (operand distribution: one tile from the temporal buffer to EVERY active PE).
//     The fabric forwards one copy along a tree, and every one of the h*w PEs receives it over
//     its own incoming link -- so h*w links carry a copy of each transaction. The previous
//     model charged the AVERAGE Manhattan hop count instead (15 on a 16x16 array), which is the
//     cost of one average unicast, not of a fanout to 256 endpoints. It understated a 16x16
//     array's distribution energy by roughly 17x. This was invisible because every shipped
//     config sets noc_energy = 0.
//   UNICAST (output write-back: each PE sends its OWN partial sum to the array edge).
//     Here each transaction really does travel one route, so the average Manhattan distance per
//     transaction IS the right multiplier -- the old formula was correct for this direction.
//
// Treating both directions with one formula is what made the distribution case wrong, so the
// mode is an explicit parameter rather than a default.
struct spatial_noc_cost_t {
    double latency_fill_hops;   // max(1, max Manhattan hops) - 1; 0 for single-hop fabrics
    double link_traversals;     // fabric links crossed per delivered transaction; see below
};

// RE6: what ONE link traversal IS, for the spatial (intra-array) NoC.
//
// `noc_energy` prices ONE traversal of ONE ROUTER-TO-ROUTER LINK INSIDE the array fabric. It does
// NOT price:
//   * the GLB <-> PE-array attach link. That link is charged separately by the GLB section's
//     `transfer_energy` (the report's "PEs - Global buffer" axis), so counting it here would bill
//     the same wire twice.
//   * endpoint injection/ejection into a PE. Per the component boundary the model documents -- a
//     transfer charges the source's read to the SOURCE, the link to its OWNER, and the
//     destination's write to the DESTINATION -- that is the PE local buffer's write energy.
//
// Under that contract the counts are edge sets, not receiver counts:
//   MULTICAST (operand distribution, one tile to every active PE): a spanning tree over the N
//     active routers, rooted at the router the GLB attaches to, has exactly N-1 edges. So the
//     answer to the "N or N-1" question is N-1. N would be right only under a contract that also
//     prices the attach link, which this model prices elsewhere.
//   UNICAST (partial-sum write-back, each PE to the array edge): the route from a PE's router to
//     the corner router, averaged over the active PEs -- (h-1)/2 + (w-1)/2 internal hops.
// A 1x1 mesh therefore has ZERO internal links: its data arrives entirely over the attach link.
// A BUS or CROSSBAR is one shared medium: one traversal, no route depth.
//
// spatial_noc_link_contract() returns this contract as one line for the report, so a reader can
// tell which convention produced the number (RE6 completion condition 4).
const char *spatial_noc_link_contract(noc_type_t topology, bool multicast);

// RE6 condition 3: build the ACTUAL multicast edge set for an active mesh shape and return its
// size, by enumerating the dimension-order spanning tree rather than evaluating a closed form.
// spatial_noc_cost() must agree with it for every shape, which is what stops the closed form from
// drifting away from a physical link count.
unsigned spatial_multicast_edge_count(unsigned active_height, unsigned active_width);

// RE6: the same statement for the NoP. Its contract DIFFERS from the array's on purpose: the NoP's
// source is the multi-chip staging buffer, whose link into chip 0 is part of the NoP itself and is
// priced by no other axis -- so that ingress link IS counted here, and nop_delivery_cost() adds it
// explicitly. Each chip's own delivery is its GLB's write energy, as always.
const char *nop_link_contract(noc_type_t topology, bool multicast);

spatial_noc_cost_t spatial_noc_cost(noc_type_t topology, unsigned active_height,
                                    unsigned active_width, bool multicast);

// SY2/L9: fill and DRAIN of a 2D systolic array of the given ACTIVE shape.
//
// `skew_fill_hops` is the operand skew -- the farthest PE starts (h-1)+(w-1) hop cycles
// after the first, because operands ripple through one pipeline register per hop. This is
// the injection wavefront the load path already charges via spatial_noc_cost(MESH, ...).
//
// `drain_hops` is what an injection-only model misses: after the LAST operand enters, the
// last partial sum still has to propagate along the accumulation dimension to reach the
// edge accumulator, and until it has, the array cannot switch to a different weight
// residency. That makes it a tile-boundary bubble, not a one-time layer cost. Partial sums
// travel down the columns, so the depth is the active HEIGHT, not the diameter.
// (Gemmini's RTL-measured ~DIM cycles of load/accumulate-drain bubble per array-wide weight
// residency is this quantity, which is why a config that calibrates
// weight_fold_fill_cycle uses the measurement instead of this estimate.)
struct systolic_pipeline_cost_t {
    double skew_fill_hops;
    double drain_hops;
};

systolic_pipeline_cost_t systolic_pipeline_cost(unsigned active_height, unsigned active_width);

// Multi-chip NoP: a serialized bus fabric, or a routed-unicast mesh where each chip's
// stream traverses its Manhattan distance from the package ingress at grid (0,0).
bool is_supported_multi_chip_nop(noc_type_t topology);

// Manhattan hop count from the package ingress (0,0) to chip `chip_index`, laid out
// row-major on a grid with `grid_width` columns. At least one link is traversed.
unsigned nop_route_hops(unsigned chip_index, unsigned grid_width);

// L7: cost of delivering ONE tile over the NoP to `active_chips` chips laid out row-major on
// a grid `grid_width` wide, entering the package at a single ingress next to chip (0,0).
//
// The ingress is the shared bottleneck: every delivery crosses it. That is what separates
// the two delivery contracts.
//   * ROUTED UNICAST -- each chip needs its own distinct chunk, so the ingress carries one
//     copy per chip and those copies serialize there.
//   * MULTICAST -- every chip needs the SAME tile (a broadcast datatype: one whose tile does
//     not depend on the CHIPS_X/CHIPS_Y-split dimension), so the fabric forwards ONE copy
//     that fans out along a tree. The ingress carries one copy no matter how many chips
//     receive it. On a bus this is not even an optimization: one transmission is physically
//     seen by every receiver, and charging it per chip was simply wrong.
//
// `bottleneck_link_tiles` is how many copies cross the busiest link (the ingress) and so
// sets the fabric's serialized busy time; `total_link_traversals` is the sum over all links
// and so sets energy (every hop transfer burns link energy); `fill_hops` is the route depth
// beyond the first hop, charged ONCE per delivery because the hops of concurrent routes
// pipeline rather than queueing behind each other.
struct nop_delivery_cost_t {
    double bottleneck_link_tiles;
    double total_link_traversals;
    double fill_hops;
};

nop_delivery_cost_t nop_delivery_cost(noc_type_t topology, unsigned grid_width,
                                      unsigned active_chips, bool multicast);

// B12: derive a link bit-width (bits) from bandwidth (GB/s) and frequency (MHz),
// warning on stderr when the result truncates fractionally -- silent truncation
// biases every ceil(bits/bitwidth) transaction count. Returns 0 when frequency <= 0.
unsigned derived_link_bitwidth(const char *component, double bandwidth, double frequency);

// Balanced 2-input adder-tree reduction over `leaves` active operands.
// depth is the number of sequential levels (parallel additions per level);
// num_additions is the total adder-operation count across the whole tree.
struct adder_tree_reduction_cost_t {
    unsigned num_additions;
    unsigned depth;
};

adder_tree_reduction_cost_t adder_tree_reduction_cost(unsigned leaves);

#endif
