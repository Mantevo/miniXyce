
import random
from collections import namedtuple


class Node:
    
    def __init__(self, name):
        self.name = name


class PI_Node(Node):

    def __init__(self, name):
        super().__init__(name)


Op = namedtuple('Op', ['name', 'aryness'])


class Intermediate_Node(Node):

    def __init__(self, name, op, args):
        super().__init__(name)
        self.op = op
        self.args = args


class PO_Node(Node):

    def __init__(self, name, arg):
        super().__init__(name)
        self.arg = arg


class Random_Function:

    def __init__(self, in_vec_len, out_vec_len, num_ops):
        self._num_PIs = in_vec_len
        self._num_POs = out_vec_len
        self._num_intermediate_nodes = num_ops
        self._prefix = 'z'
        self._generate_graph()

    def _generate_graph(self):
        num_nodes = 0
        self._PIs = []
        for idx in range(self._num_PIs):
            node_name = self._get_node_name(num_nodes)
            PI = PI_Node(node_name)
            self._PIs.append(PI)
            num_nodes += 1
        unused_nodes = set(PI for PI in self._PIs)
        self._intermediate_nodes = []
        for idx in range(self._num_intermediate_nodes):
            node_name = self._get_node_name(num_nodes)
            u = self._create_random_intermediate_node(node_name, unused_nodes)
            self._intermediate_nodes.append(u)
            num_nodes += 1
        self._POs = []
        for idx in range(self._num_POs):
            node_name = self._get_node_name(num_nodes)
            PO = self._create_random_PO_node(node_name, unused_nodes)
            self._POs.append(PO)
            num_nodes += 1
    
    def _get_node_name(self, node_idx):
        return '{}{}'.format(self._prefix, node_idx)

    def _create_random_intermediate_node(self, node_name, unused_nodes):
        op = self._get_random_op()
        args = []
        for idx in range(op.aryness):
            if unused_nodes:
                arg = random.choice(list(unused_nodes))
                unused_nodes.remove(arg)
            else:
                arg = random.choice(self._PIs + self._intermediate_nodes)
            args.append(arg)
        ans = Intermediate_Node(node_name, op, args)
        unused_nodes.add(ans)
        return ans

    def _get_random_op(self, safe=True):
        s = '_safe' if safe else ''
        return random.choice(getattr(self, '_get{}_ops'.format(s))())

    def _get_ops(self):
        ans = [ Op('UNARY_PLUS', 1), 
                Op('UNARY_MINUS', 1), 
                Op('PLUS', 2), 
                Op('MINUS', 2), 
                Op('STAR', 2), 
                Op('SLASH', 2), 
                Op('POW', 2), 
                Op('SQUARE', 1), 
                Op('CUBE', 1), 
                Op('SQRT', 1), 
                Op('CBRT', 1), 
                Op('EXP', 1), 
                Op('LOG', 1), 
                Op('LOG10', 1), 
                Op('LOG2', 1), 
                Op('ABS', 1), 
                Op('SIN', 1), 
                Op('COS', 1), 
                Op('TAN', 1), 
                Op('ASIN', 1), 
                Op('ACOS', 1), 
                Op('ATAN', 1), 
                Op('SINH', 1), 
                Op('COSH', 1), 
                Op('TANH', 1), 
                Op('ASINH', 1), 
                Op('ACOSH', 1), 
                Op('ATANH', 1), ]
        return ans

    def _get_safe_ops(self):
        ans = [ Op('UNARY_PLUS', 1), 
                Op('UNARY_MINUS', 1), 
                Op('PLUS', 2), 
                Op('MINUS', 2), 
                Op('STAR', 2), 
                Op('SQUARE', 1), 
                Op('CUBE', 1), 
                Op('CBRT', 1), 
                Op('ABS', 1), 
                Op('SIN', 1), 
                Op('COS', 1), 
                Op('ATAN', 1), 
                Op('TANH', 1), ]
        return ans

    def _create_random_PO_node(self, node_name, unused_nodes):
        if unused_nodes:
            arg = random.choice(list(unused_nodes))
            unused_nodes.remove(arg)
        else:
            arg = random.choice(self._PIs + self._intermediate_nodes)
        ans = PO_Node(node_name, arg)
        return ans

    def export_code(self, func_name, in_vec_name, out_vec_name, language):
        assert language in ['Python', 'C++'], 'ERROR: Unrecognized language'
        export_func_name = { 'Python': '_export_python_code', 
                             'C++': '_export_cpp_code' }[language]
        export_func = getattr(self, export_func_name)
        export_func(func_name, in_vec_name, out_vec_name)

    def _export_python_code(self, func_name, in_vec_name, out_vec_name):

        print('import math')

        print('')
        print('def {}({}):'.format(func_name, in_vec_name))

        print('')
        for idx, PI in enumerate(self._PIs):
            print('    {} = {}[{}]'.format(PI.name, in_vec_name, idx))

        print('')
        for u in self._intermediate_nodes:
            line = '    {} = '.format(u.name)
            if u.op.name == 'UNARY_PLUS':
                line += '+{}'.format(u.args[0].name)
            elif u.op.name == 'UNARY_MINUS':
                line += '-{}'.format(u.args[0].name)
            elif u.op.name == 'PLUS':
                line += '{} + {}'.format(u.args[0].name, u.args[1].name)
            elif u.op.name == 'MINUS':
                line += '{} - {}'.format(u.args[0].name, u.args[1].name)
            elif u.op.name == 'STAR':
                line += '{} * {}'.format(u.args[0].name, u.args[1].name)
            elif u.op.name == 'SLASH':
                line += '{} / {}'.format(u.args[0].name, u.args[1].name)
            elif u.op.name == 'POW':
                line += '{} ** {}'.format(u.args[0].name, u.args[1].name)
            elif u.op.name == 'SQUARE':
                line += '{} ** 2.0'.format(u.args[0].name)
            elif u.op.name == 'CUBE':
                line += '{} ** 3.0'.format(u.args[0].name)
            elif u.op.name == 'CBRT':
                line += '{} ** (1.0/3.0)'.format(u.args[0].name)
            elif u.op.name == 'ABS':
                line += 'std::abs({})'.format(u.args[0].name)
            elif u.op.name in ( 'SQRT', 'EXP', 'LOG', 'LOG10', 'LOG2', 'SIN', 
                'COS', 'TAN', 'ASIN', 'ACOS', 'ATAN', 'SINH', 'COSH', 'TANH',
                'ASINH', 'ACOSH', 'ATANH' ):
                line += 'math.{}({})'.format(u.op.name.lower(), u.args[0].name)
            else:
                assert False, 'ERROR: Unrecognized operation'
            print(line)

        print('')
        for PO in self._POs:
            print('    {} = {}'.format(PO.name, PO.arg.name))

        print('')
        print('    {} = []'.format(out_vec_name))
        for PO in self._POs:
            print('    {}.append({})'.format(out_vec_name, PO.name))

        print('')
        print('    return {}'.format(out_vec_name))
    
    def _export_cpp_code(self, func_name, in_vec_name, out_vec_name):

        print('#include <cmath>')
        print('#include <vector>')

        print('')
        print('std::vector<double> {}(std::vector<double> &{})'.format(
            func_name, in_vec_name))

        print('{')

        for idx, PI in enumerate(self._PIs):
            print('    double {} = {}[{}];'.format(PI.name, in_vec_name, idx))

        print('')
        for u in self._intermediate_nodes:
            line = '    double {} = '.format(u.name)
            if u.op.name == 'UNARY_PLUS':
                line += '+{}'.format(u.args[0].name)
            elif u.op.name == 'UNARY_MINUS':
                line += '-{}'.format(u.args[0].name)
            elif u.op.name == 'PLUS':
                line += '{} + {}'.format(u.args[0].name, u.args[1].name)
            elif u.op.name == 'MINUS':
                line += '{} - {}'.format(u.args[0].name, u.args[1].name)
            elif u.op.name == 'STAR':
                line += '{} * {}'.format(u.args[0].name, u.args[1].name)
            elif u.op.name == 'SLASH':
                line += '{} / {}'.format(u.args[0].name, u.args[1].name)
            elif u.op.name == 'POW':
                line += 'std::pow({}, {})'.format(u.args[0].name, 
                                                  u.args[1].name)
            elif u.op.name == 'SQUARE':
                line += 'std::pow({}, 2.0)'.format(u.args[0].name)
            elif u.op.name == 'CUBE':
                line += 'std::pow({}, 3.0)'.format(u.args[0].name)
            elif u.op.name == 'ABS':
                line += 'std::abs({})'.format(u.args[0].name)
            elif u.op.name in ( 'SQRT', 'CBRT', 'EXP', 'LOG', 'LOG10', 'LOG2', 
                'SIN', 'COS', 'TAN', 'ASIN', 'ACOS', 'ATAN', 'SINH', 'COSH', 
                'TANH', 'ASINH', 'ACOSH', 'ATANH' ):
                line += 'std::{}({})'.format(u.op.name.lower(), u.args[0].name)
            else:
                assert False, 'ERROR: Unrecognized operation'
            line += ';'
            print(line)

        print('')
        for PO in self._POs:
            print('    double {} = {};'.format(PO.name, PO.arg.name))

        print('')
        print('    std::vector<double> {};'.format(out_vec_name))
        for PO in self._POs:
            print('    {}.push_back({});'.format(out_vec_name, PO.name))

        print('')
        print('    return {};'.format(out_vec_name))

        print('}')
    


if __name__ == '__main__':
    
    # define parameters
    func_name = 'diode_l1'
    in_vec_name, in_vec_len = 'x', 3
    out_vec_name, out_vec_len = 'y', 3
    num_ops = 50
    language = 'C++' # one of 'C++', 'Python'

    # generate a random function graph
    func = Random_Function(in_vec_len, out_vec_len, num_ops)
    func.export_code(func_name, in_vec_name, out_vec_name, language)

