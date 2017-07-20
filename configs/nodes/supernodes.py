from astrojc.pipeline import *

__all__ = ['ForSuperNode']

class ForSuperNode(Node):
    '''This is a temporary node that will act as node group with for function.'''
    '''This time, the nodes have to be created outside this.'''
    def __init__(self, name, subcontext = None, **kwargs):
        Node.__init__(self, name, **kwargs)
        if isinstance(subcontext, Context):
            self._subctx = subcontext
        elif subcontext == None:
            self._subctx = Context(name)
        else:
            raise ValueError('Invalid subcontext.')

        self._setup_subctx()

    def _setup_subctx(self):
        self._subctx.log = self._ctx.log
        self.add_node = self._subctx.add_node
        self.del_node = self._subctx.del_node
        self.create_node = self._subctx.create_node

    def setup(self):
        self.add_input('[i]')
        self.add_output('i')
        self.add_input('[j]', optional = True)
        self.add_output('j')
        self.add_input('[k]', optional = True)
        self.add_output('k')

        self.add_input('result', optional = True)
        self.add_output('[results]')

    def run(self):
        iter_size = len(self['input.[i]'])

        have_j = False
        have_k = False

        if self.get_dock('input.[j]').have_value:
            if len(self['input.[j]']) != iter_size:
                raise ValueError('[i] and [j] list do not match in size.')
            else:
                have_j = True

        if self.get_dock('input.[k]').have_value:
            if len(self['input.[k]']) != iter_size:
                raise ValueError('[i] and [k] list do not match in size.')
            else:
                have_k = True

        results = [None] * iter_size

        for i in range(iter_size):
            self.log.debug('Running {} {}/{}'.format(self.name, i+1, iter_size))
            self['output.i'] = self['input.[i]'][i]
            if have_j:
                self['output.j'] = self['input.[j]'][i]
            if have_k:
                self['output.k'] = self['input.[k]'][i]

            self._subctx.run()

            results[i] = self['input.result']
            self['output.i'] = empty_value
            self['output.j'] = empty_value
            self['output.k'] = empty_value
            for j in self._subctx._nodes.values():
                j.set_state('idle')

        self['output.[results]'] = results
        for i in self._subctx._nodes.values():
            i.set_state('done')

    def _check_ready(self):
        self['output.i'] = self['input.[i]']
        self['output.j'] = self['input.[j]']
        self['output.k'] = self['input.[k]']
        return Node._check_ready(self)
        ## I think this is a source of bug and useless
        #for i in self._nodes.values():
        #    if not i.is_ready:
        #        return False

    def get_dock(self, key):
        '''Gets the object of a dock with a given key.'''
        try:
            return Node.get_dock(self, key)
        except:
            try:
                return self._subctx.get_dock(key)
            except:
                raise ValueError('{} key doesn\'t correspond to any '
                                 'dock in this supernode.'.format(key))

    def items(self):
        items = Node.vars(self)

        for i,v in self._subctx.vars():
            items[i] = v

        return items

    def __getitem__(self, key):
        try:
            return Node.__getitem__(self, key)
        except:
            try:
                return self._subctx[key]
            except:
                raise ValueError('{} key doesn\'t correspond to any property '
                                 'or node from this supernode.'.format(key))

    def __setitem__(self, key, value):
        try:
            Node.__setitem__(self, key, value)
        except:
            try:
                self._subctx[key] = value
            except:
                raise ValueError('{} key doesn\'t correspond to any property '
                                 'or node from this supernode.'.format(key))
