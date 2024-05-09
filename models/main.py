import os
import sys
pwd_path = os.path.dirname(os.path.abspath(os.path.dirname(__file__)))
util_path = os.path.join(pwd_path, 'python_util')
sys.path.append(util_path)
import argparse
import network_builder


def network_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('-n', 
                        metavar='network.cfg', 
                        type=str,
                        help='Network configuration file', default='NONE',
                        )
    return parser

def build_network(m_network):
    t_network = m_network.replace('networks/', '')
    #print(t_network)
    parser = network_parser()
    args = parser.parse_args()
    network = t_network

    DNN_model = network_builder.network_builder(network)


if __name__ == '__main__':

    print('run Python main function')
    parser = network_parser()
    args = parser.parse_args()
    network = args.n

    DNN_model = network_builder.network_builder(network)
