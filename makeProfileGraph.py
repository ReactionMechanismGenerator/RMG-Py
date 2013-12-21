#!/usr/bin/env python
# -*- coding: utf-8 -*-

def makeProfileGraph(stats_file,thresh_node,thresh_edge):
    """
    Uses gprof2dot to create a graphviz dot file of the profiling information.
    
    This requires the gprof2dot package available via `pip install gprof2dot`.
    Renders the result using the program 'dot' via a command like
    `dot -Tpdf input.dot -o output.pdf`.
    """
    try:
        from external.gprof2dot import gprof2dot
    except ImportError:
        try:
            from external import gprof2dot
        except ImportError:    
            print('Package gprof2dot not found. Unable to create a graph of the profile statistics.')
            print("`pip install gprof2dot` if you don't have it.")
            return
    import subprocess
    m = gprof2dot.Main()
    class Options:
        pass
    m.options = Options()
    m.options.node_thres = thresh_node# default 0.8
    m.options.edge_thres = thresh_edge # default 0.1
    m.options.strip = False
    m.options.wrap = True
    m.theme = m.themes['color'] # bw color gray pink
    parser = gprof2dot.PstatsParser(stats_file)
    m.profile = parser.parse()
    dot_file = stats_file + '.dot'
    m.output = open(dot_file,'wt')
    m.write_graph()
    m.output.close()
    try:
        subprocess.check_call(['dot', '-Tpdf', dot_file, '-o', '{0}.pdf'.format(dot_file)])
    except subprocess.CalledProcessError:
        print("Error returned by 'dot' when generating graph of the profile statistics.")
        print("To try it yourself:\n     dot -Tpdf {0} -o {0}.pdf".format(dot_file))
    except OSError:
        print("Couldn't run 'dot' to create graph of profile statistics. Check graphviz is installed properly and on your path.")
        print("Once you've got it, try:\n     dot -Tpdf {0} -o {0}.pdf".format(dot_file))
    else:
        print("Graph of profile statistics saved to: \n {0}.pdf".format(dot_file))
        
if __name__ == '__main__':
 
    import argparse
     
    parser = argparse.ArgumentParser(description="Creates a call graph with profiling information.")
    parser.add_argument('FILE', type=str, default='RMG.profile',nargs='?', help='.profile file (default file is RMG.profile)')
    parser.add_argument('THRESH_NODE', type=float, default=0.8,nargs='?', help='threshold percentage value for nodes (default value is 0.8)')
    parser.add_argument('THRESH_EDGE', type=float, default=0.1, nargs='?', help='threshold percentage value for nodes (default value is 0.1)') 
    args = parser.parse_args()
    stats_file=args.FILE
    thresh_node=args.THRESH_NODE
    thresh_edge=args.THRESH_EDGE
    
    makeProfileGraph(stats_file,thresh_node,thresh_edge)
    
    
