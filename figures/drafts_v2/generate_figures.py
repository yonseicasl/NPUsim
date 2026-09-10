#!/usr/bin/env python3
"""Editable, monochrome paper figures; no external assets or network access."""
from pathlib import Path
from html import escape
import xml.etree.ElementTree as ET

OUT = Path(__file__).resolve().parent
FONT = 'Arial, Helvetica, Liberation Sans, sans-serif'
class SVG:
    def __init__(self,w,h,title,desc):
        self.w,self.h=w,h
        self.parts=[f'<svg xmlns="http://www.w3.org/2000/svg" width="{w}" height="{h}" viewBox="0 0 {w} {h}" role="img" aria-labelledby="title desc">',f'<title id="title">{escape(title)}</title><desc id="desc">{escape(desc)}</desc>', '''<defs>
<marker id="arrow" viewBox="0 0 8 8" refX="7" refY="4" markerWidth="6" markerHeight="6" orient="auto-start-reverse"><path d="M0 0 L8 4 L0 8 Z" fill="#111"/></marker>
<pattern id="stall" width="7" height="7" patternUnits="userSpaceOnUse"><rect width="7" height="7" fill="white"/><path d="M-2 2 L2 -2 M0 7 L7 0 M5 9 L9 5" stroke="#777" stroke-width="1"/></pattern>
</defs>''',f'<g font-family="{FONT}" font-size="18" fill="#111" stroke-linejoin="miter">',f'<rect width="{w}" height="{h}" fill="white"/>']
    def raw(self,x): self.parts.append(x)
    def rect(self,x,y,w,h,fill='white',dash=False,sw=1.3):
        self.raw(f'<rect x="{x}" y="{y}" width="{w}" height="{h}" fill="{fill}" stroke="#111" stroke-width="{sw}"'+(' stroke-dasharray="6 4"' if dash else '')+'/>')
    def text(self,x,y,s,size=18,bold=False,anchor='start',italic=False,fill='#111',mono=False):
        self.raw(f'<text x="{x}" y="{y}" font-size="{size}" text-anchor="{anchor}"'+(' font-weight="bold"' if bold else '')+(' font-style="italic"' if italic else '')+(' font-family="Liberation Mono, monospace"' if mono else '')+f' fill="{fill}">{escape(s)}</text>')
    def lines(self,x,y,lines,size=18,dy=24,**kw):
        for i,s in enumerate(lines):self.text(x,y+i*dy,s,size,**kw)
    def path(self,pts,arrow=False,dash=False,sw=1.5,both=False):
        d='M'+' L'.join(f'{x},{y}' for x,y in pts)
        self.raw(f'<path d="{d}" fill="none" stroke="#111" stroke-width="{sw}"'+(' stroke-dasharray="6 4"' if dash else '')+(' marker-end="url(#arrow)"' if arrow else '')+(' marker-start="url(#arrow)"' if both else '')+'/>')
    def box(self,x,y,w,h,title,lines=(),size=18,fill='#ededed',dash=False):
        self.rect(x,y,w,h,fill,dash)
        self.text(x+12,y+25,title,size,bold=True)
        self.lines(x+12,y+49,lines,size-1,dy=23)
    def group(self,x,y,w,h,title,size=20):
        self.rect(x,y,w,h,'none',True)
        tw=len(title)*size*.54+16
        self.raw(f'<rect x="{x+10}" y="{y-12}" width="{tw}" height="22" fill="white"/>')
        self.text(x+16,y+5,title,size,bold=True)
    def pe(self,x,y,w=27,h=27,label='PE',fill='white',size=14):
        self.rect(x,y,w,h,fill,sw=1)
        self.text(x+w/2,y+h/2+size*.34,label,size,anchor='middle')
    def finish(self,name):
        s='\n'.join(self.parts+['</g></svg>'])
        ET.fromstring(s)
        (OUT/name).write_text(s)
        return s

# Fig. 2: input / main / output, with separate hardware and scheduling inputs.
s=SVG(1220,500,'Fig. 2. NPUsim overview','Separate hardware specifications, scheduling schemes and DNN workloads drive decoupled execution schedulers and modular hardware models. Component activities determine timing, utilization and energy; functional execution is optional.')
s.group(12,24,277,459,'Input module')
s.group(344,24,607,459,'Main module')
s.group(999,24,209,459,'Output module')
s.box(26,55,248,100,'Hardware specification',['PEs, buffers, links, DRAM','Structure, precision, format'])
s.box(26,188,248,79,'Scheduling scheme',['Dataflow & mapping'])
s.box(26,310,248,105,'DNN workload / PyTorch',['Layer dimensions','Tensor values (optional)'])
s.box(26,434,248,35,'Accelergy / CACTI',size=18,fill='white')
s.box(365,188,230,128,'Execution schedulers',['Per-boundary tile schedules','Tile offsets & reuse counts','Transfer & compute order'])
s.group(684,52,250,347,'Simulation engine',19)
s.box(698,81,222,83,'Component models',['Instantiate & connect','PE / buffer / link / DRAM'],size=18)
s.box(698,196,222,98,'Cycle-level timing',['Dependencies / contention','Buffering / overlap'],size=18)
s.box(698,327,222,54,'Functional execution',['Actual tensor arithmetic'],size=17,fill='white',dash=True)
s.box(365,431,230,37,'Component cost model',size=18,fill='white')
s.box(698,431,222,37,'Activity / cost accounting',size=17,fill='#ddd')
s.path([(274,107),(698,107)],True)
s.path([(274,227),(365,227)],True)
s.path([(274,346),(319,346),(319,284),(365,284)],True)
s.text(295,302,'Layers',16)
s.path([(274,388),(609,388)],False,True)
s.path([(621,388),(657,388),(657,355),(698,355)],True,True)
s.text(450,379,'Tensor values (optional)',17,anchor='middle')
s.path([(595,217),(698,217)],True)
s.text(646,207,'Tiles / order',16,anchor='middle')
s.path([(698,278),(595,278)],True,True)
s.text(646,266,'Requests',16,anchor='middle')
s.path([(809,164),(809,196)],True)
s.path([(809,294),(809,327)],True,True)
s.path([(920,246),(941,246),(941,411),(809,411),(809,431)],True)
s.path([(809,381),(809,411)],False,True)
s.path([(274,451),(365,451)],True)
s.path([(595,448),(698,448)],True)
s.text(646,436,'Unit costs',16,anchor='middle')
s.path([(615,448),(615,306),(698,306),(698,294)],True)
s.path([(920,450),(973,450),(973,223),(1013,223)],True)
s.box(1013,64,181,304,'Runtime statistics',size=18,fill='white')
s.lines(1025,119,['Execution cycles','Energy','PE utilization','Buffer occupancy','Operation count','OPS / intensity'],18,dy=36)
s.box(1013,391,181,77,'DNN output',['Inference accuracy'],size=18,fill='white',dash=True)
s.text(1103,458,'(functional mode)',16,anchor='middle',italic=True)
s.path([(973,429),(1013,429)],True,True)
s.text(365,345,'What moves, and in what order',17,italic=True)

s.finish('fig2_npusim_overview.svg')

# Fig. 3: declarative composition with compact, distinct architecture examples.
s=SVG(620,797,'Fig. 3. Configurable architecture modeling','Hardware configuration blocks select and parameterize a shared component library. The same library composes spatial, systolic and multi-chip designs, and supports optional optimization modules.')
s.box(12,24,209,235,'Hardware specification',size=17.5,fill='white')
s.lines(25,80,['Component blocks','• Type / organization','• Dimensions / capacity','• Bandwidth / clock','• Access cycles / energy','• Precision / format'],17,dy=29)
s.path([(221,140),(262,140)],True)
s.text(241,119,'select',15,anchor='middle',italic=True)
s.group(266,24,342,394,'Component library',19)
rows=[(64,'PE','MAC lanes + I / W / O local buffers'),(127,'Array / interconnect','Spatial / systolic / tree; bus / mesh'),(190,'Global buffer','Shared / partitioned; double / bypass'),(253,'Package network','Chip grid + multicast delivery'),(316,'DRAM','Analytical model / DRAMsim3')]
for y,title,detail in rows:
 s.text(281,y,title,18,bold=True)
 s.text(281,y+25,detail,16.5)
 if y<316:s.path([(280,y+43),(593,y+43)],sw=.7)
s.rect(279,365,315,39,'#ededed',True)
s.text(436,390,'Optimization modules',17,anchor='middle',bold=True)
s.box(12,290,209,128,'Optional module blocks',['Weight decompressor','SFU / KV-cache unit'],size=17,fill='white',dash=True)
s.path([(221,385),(279,385)],True,True)
s.path([(437,418),(437,429),(310,429)],sw=1.3)
s.text(310,450,'Instantiate & connect at initialization',18,bold=True,anchor='middle')
s.path([(310,456),(310,468),(102,468),(102,494)],True)
s.path([(310,468),(310,494)],True)
s.path([(310,468),(518,468),(518,494)],True)
# Archetypes use the same primitives, with no numeric claims about real chips.
for x,title,sub in [(12,'(a) Spatial','Eyeriss-like'),(220,'(b) Systolic','TPU-like'),(428,'(c) Multi-chip','Simba-like')]:
 s.rect(x,494,180,278,'white')
 s.text(x+90,519,title,18,bold=True,anchor='middle')
 s.text(x+90,542,sub,17,italic=True,anchor='middle')
 s.rect(x+17,558,146,28,'#ddd')
 s.text(x+90,578,'HBM' if x==220 else 'DRAM',17,anchor='middle')
 s.path([(x+90,586),(x+90,606)],True,both=True)
# Spatial: per tensor global buffers / bus / PE array.
x=12
for i,l in enumerate(['I','W','O']):
 s.rect(x+17+i*49,607,48,28,'#eee'); s.text(x+41+i*49,627,l,16,anchor='middle')
s.text(x+90,654,'Partitioned buffer',16,anchor='middle')
s.path([(x+90,657),(x+90,671)])
s.path([(x+24,671),(x+156,671)],sw=2)
for j in range(3):
 px=x+28+49*j
 s.path([(px+13,671),(px+13,690)],True)
 s.pe(px,690,27,27)
 s.path([(px+13,717),(px+13,731)],True)
 s.pe(px,731,27,27)
s.text(x+174,687,'Bus',14,anchor='end')
# Systolic: shared buffer and 3 x 3 mesh.
x=220
s.rect(x+17,607,146,28,'#eee');s.text(x+90,627,'Shared buffer',16,anchor='middle')
for r in range(3):
 for c in range(3):
  px=x+31+c*47;py=654+r*36
  if r==0:s.path([(px+13,635),(px+13,py)],True)
  if c<2:s.path([(px+27,py+13),(px+47,py+13)],True)
  if r<2:s.path([(px+13,py+27),(px+13,py+36)],True)
  s.pe(px,py,27,27)
# Multi-chip: a 2 x 2 package mesh with buffer and PE structure in each chip.
x=428
s.rect(x+12,608,156,149,'white',True)
s.text(x+90,630,'Package network',16,anchor='middle')
for r in range(2):
 for c in range(2):
  px=x+22+c*77;py=643+r*59
  s.rect(px,py,59,44,'#eee');s.text(px+29,py+18,'Chip',16,anchor='middle')
  s.rect(px+7,py+25,15,12,'white',sw=.7)
  for q in range(3):s.rect(px+27+8*q,py+25,6,12,'white',sw=.7)
  if c==0:s.path([(px+59,py+22),(px+77,py+22)],True,both=True)
  if r==0:s.path([(px+29,py+44),(px+29,py+59)],True,both=True)
s.text(310,791,'I / W / O: input / weight / output     SFU: special function unit',15.5,anchor='middle')
s.finish('fig3_configurable_architecture.svg')

# Fig. 4: two panels, also exported individually for LaTeX subfloat labels.
def schedule_panel(s):
 s.text(12,24,'(a) Decoupled execution scheduling',21,bold=True)
 s.box(12,43,273,64,'Scheduling scheme',['Dataflow + mapping at each level'],size=18,fill='white')
 s.box(318,43,270,64,'DNN layer dimensions',size=18,fill='white')
 s.path([(149,107),(149,125),(300,125),(300,143)],True)
 s.path([(453,107),(453,125),(300,125)])
 s.box(12,143,576,151,'Per-boundary tile schedules',size=19,fill='#eee')
 for y in [181,208,235,262]:s.path([(24,y),(576,y)],sw=.65)
 for x in [143,395]:s.path([(x,181),(x,285)],sw=.65)
 s.text(32,201,'Tensor',18,bold=True);s.text(159,201,'Tile offsets (ordered)',18,bold=True);s.text(409,201,'Reuse counts',18,bold=True)
 s.text(32,228,'Weight',18);s.text(159,228,'w₀, w₁, …',18);s.text(409,228,'4, 4, …',18)
 s.text(32,255,'Input',18);s.text(159,255,'i₀, i₁, …',18);s.text(409,255,'1, 1, …',18)
 s.text(32,282,'Output',18);s.text(159,282,'o₀, o₁, …',18);s.text(409,282,'1, 1, …',18)
 s.text(12,324,'WS example at the local buffer',19,bold=True)
 s.text(12,346,'One weight tile is reused across four input tiles.',17)
 s.text(16,379,'Weight',18);s.text(16,417,'Input',18)
 for i in range(4):
  s.rect(121+i*90,393,90,33,'white');s.text(166+i*90,416,f'i{chr(0x2080+i)}',18,anchor='middle')
 s.rect(121,357,360,30,'#ccc');s.text(301,379,'w₀ held in place (reuse = 4)',18,anchor='middle')
 s.rect(504,357,70,30,'white');s.text(539,379,'w₁',18,anchor='middle')
 s.path([(481,372),(504,372)],True)
 s.text(539,416,'…',22,anchor='middle')
 s.path([(121,440),(574,440)],True);s.text(346,462,'Logical order (not elapsed cycles)',16,anchor='middle',italic=True)
 s.text(12,497,'Demand-driven execution on hardware',19,bold=True)
 # Hardware nodes with boundary schedulers below each link.
 nodes=[(12,96,'DRAM'),(132,104,'Package'),(260,86,'GB'),(370,86,'LB'),(480,108,'MAC / RF')]
 for x,w,l in nodes:
  s.rect(x,555,w,38,'#eee');s.text(x+w/2,580,l,18,anchor='middle')
 for i in range(4):
  a=nodes[i][0]+nodes[i][1];b=nodes[i+1][0]
  s.path([(a,574),(b,574)],True)
  mid=(a+b)/2
  s.rect(mid-19,607,38,27,'white');s.text(mid,627,f'S{chr(0x2080+i)}',17,anchor='middle')
  s.path([(mid,607),(mid,585)],True)
 for i in range(4):
  upstream=nodes[i][0]+nodes[i][1]/2;downstream=nodes[i+1][0]+nodes[i+1][1]/2
  s.path([(downstream,552),(downstream,533),(upstream,533),(upstream,552)],True,True)
 s.text(299,525,'Request escalates if the upstream level is exhausted',16,anchor='middle')
 s.text(300,657,'S₀–S₃: scheduler at each memory boundary; solid arrows: delivery',16,anchor='middle')

def timing_panel(s):
 s.text(12,24,'(b) Cycle-level timing and energy',21,bold=True)
 s.box(12,43,576,64,'Hardware determines elapsed time',['Access latency, link bandwidth, resource availability'],size=19,fill='white')
 s.text(12,139,'Illustrative pipeline: two buffer slots per boundary',18,bold=True)
 # Explicit buffering in the compact stage sketch.
 for x,label,w in [(12,'DRAM',105),(242,'GB',92),(469,'LB / MACs',119)]:
  s.rect(x,160,w,38,'#eee');s.text(x+w/2,185,label,17,anchor='middle')
 for a,b in [(117,242),(334,469)]:
  s.path([(a,179),(b,179)],True)
  cx=(a+b)/2
  for k in range(2):s.rect(cx-26+k*27,167,23,24,'white')
  s.text(cx,220,'2 slots',16,anchor='middle')
 # Timing example, normalized integer cycle costs: memory 2, transfer 2, compute 4.
 x0=157; scale=20.6; y0=271
 s.path([(x0,y0),(x0+20*scale+11,y0)],True,sw=1)
 for t in range(0,21,2):
  x=x0+t*scale
  s.path([(x,y0),(x,462)],dash=True,sw=.45)
  s.text(x,258,str(t),16,anchor='middle')
 s.text(368,241,'Example cycles',16,anchor='middle')
 fills=['white','#ddd','#aaa','#707070']
 rows=[(300,'DRAM → GB',[(0,2),(2,4),(4,6),(6,8)]),(357,'GB → LB',[(2,4),(4,6),(8,10),(12,14)]),(414,'MAC compute',[(4,8),(8,12),(12,16),(16,20)])]
 for y,label,spans in rows:
  s.text(12,y+23,label,17)
  s.path([(x0,y+34),(x0+20*scale,y+34)],sw=.7)
  for i,(a,b) in enumerate(spans):
   s.rect(x0+a*scale,y,(b-a)*scale,34,fills[i],sw=1)
   s.text(x0+(a+b)*scale/2,y+23,'T'+str(i),17,anchor='middle',fill='white' if i==3 else '#111')
 for a,b in [(6,8),(10,12)]:s.rect(x0+a*scale,357,(b-a)*scale,34,'url(#stall)',sw=1)
 s.text(379,492,'Compute runs back-to-back after startup',17,anchor='middle',bold=True)
 s.path([(x0+4*scale,462),(x0+4*scale,470),(x0+20*scale,470),(x0+20*scale,462)],sw=1)
 s.rect(12,515,25,22,'url(#stall)');s.text(47,533,'Stall: downstream buffer still occupied',17)
 s.text(12,562,'T0–T3: work tiles; 2 / 2 / 4 cycles per stage.',17)
 s.text(12,585,'Slot released after its consumer finishes.',17)
 s.rect(12,605,576,58,'#eee')
 s.text(24,629,'Energy from the same traced activities',18,bold=True)
 s.text(24,652,'Access / operation counts × unit costs + leakage × time',17)

a=SVG(600,678,'Fig. 4(a). Decoupled execution scheduling','Schedulers select tile offsets and reuse counts independently of hardware timing. A weight-stationary example holds one weight tile for four input tiles. Requests propagate up the hierarchy and per-boundary schedulers select tile deliveries.')
schedule_panel(a);a.finish('fig4a_execution_scheduling.svg')
b=SVG(600,678,'Fig. 4(b). Cycle-level timing and energy','Illustrative double-buffered pipeline with four work tiles, transfer stages of two cycles and compute stages of four cycles. The global-to-local transfer stalls while local buffer slots are occupied. The compute stage is the bottleneck. These are illustrative timings, not measured results.')
timing_panel(b);b.finish('fig4b_timing_energy.svg')
s=SVG(1220,678,'Fig. 4. Decoupled scheduling and cycle-level execution','Panel a derives tile schedules and illustrates demand-driven hardware execution. Panel b composes component costs into a buffered pipeline, exposing overlap, stalls and energy activities. Timings are illustrative, not measured.')
schedule_panel(s)
s.path([(610,12),(610,665)],dash=True,sw=.7)
s.raw('<g transform="translate(620 0)">');timing_panel(s);s.raw('</g>')
s.finish('fig4_decoupled_scheduling.svg')

# Reproducible vector PDF and high-resolution bitmap previews via the system SVG renderer.
def render():
 import ctypes as C
 from ctypes.util import find_library
 cairo=C.CDLL(find_library('cairo')); rsvg=C.CDLL(find_library('rsvg-2')); gobj=C.CDLL(find_library('gobject-2.0'))
 ptr=C.c_void_p;dbl=C.c_double
 def setup(lib,name,args,result):
  f=getattr(lib,name);f.argtypes=args;f.restype=result;return f
 pdf_create=setup(cairo,'cairo_pdf_surface_create',[C.c_char_p,dbl,dbl],ptr)
 image_create=setup(cairo,'cairo_image_surface_create',[C.c_int,C.c_int,C.c_int],ptr)
 create=setup(cairo,'cairo_create',[ptr],ptr)
 scale=setup(cairo,'cairo_scale',[ptr,dbl,dbl],None)
 finish=setup(cairo,'cairo_surface_finish',[ptr],None)
 destroy=setup(cairo,'cairo_surface_destroy',[ptr],None)
 ctx_destroy=setup(cairo,'cairo_destroy',[ptr],None)
 png=setup(cairo,'cairo_surface_write_to_png',[ptr,C.c_char_p],C.c_int)
 new=setup(rsvg,'rsvg_handle_new_from_file',[C.c_char_p,ptr],ptr)
 unref=setup(gobj,'g_object_unref',[ptr],None)
 class Rect(C.Structure): _fields_=[('x',dbl),('y',dbl),('width',dbl),('height',dbl)]
 draw=setup(rsvg,'rsvg_handle_render_document',[ptr,ptr,C.POINTER(Rect),ptr],C.c_int)
 for p in sorted(OUT.glob('fig*.svg')):
  root=ET.parse(p).getroot();w=float(root.attrib['width']);h=float(root.attrib['height'])
  handle=new(str(p).encode(),None)
  paper_width=252 if p.name.startswith('fig3') else 504
  if p.name.startswith(('fig4a','fig4b')):paper_width=247.8
  factor=paper_width/w
  surface=pdf_create(str(p.with_suffix('.pdf')).encode(),paper_width,h*factor)
  ctx=create(surface);scale(ctx,factor,factor)
  viewport=Rect(0,0,w,h)
  assert draw(handle,ctx,C.byref(viewport),None)
  ctx_destroy(ctx);finish(surface);destroy(surface)
  surface=image_create(0,int(w*2),int(h*2));ctx=create(surface);scale(ctx,2,2)
  assert draw(handle,ctx,C.byref(viewport),None)
  assert png(surface,str(p.with_suffix('.png')).encode())==0
  ctx_destroy(ctx);destroy(surface);unref(handle)
  print(p.name,'→ SVG, vector PDF, PNG')
if __name__=='__main__':render()
