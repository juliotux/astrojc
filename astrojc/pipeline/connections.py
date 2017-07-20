
from .utils import check_class

class Connection(object):
    def __init__(self, input=None, output=None):
        self._input = None
        self._output = None

        try:
            self.set_input(input)
        except:
            raise ValueError('input is not valid.')

        try:
            self.set_output(output)
        except:
            raise ValueError('output is not valid.')

    @property
    def value(self):
        return self._output.get_value()

    def is_linked_to(self, dock):
        return dock == self._input or dock == self._output

    def set_input(self, input):
        if check_class(input, 'Input') or input is None:
            self._input = input
            input.link(self)
        else:
            raise ValueError('This input is not valid.')

    def set_output(self, output):
        if check_class(output, 'Output') or output is None:
            self._output = output
            output.link(self)
        else:
            raise ValueError('This output is not valid.')

    def unlink(self):
        self._input.unlink(self)
        self._output.unlink(self)
        self._input = None
        self._output = None
