def set_limits(ax, x_min=None, x_max=None, y_min=None, y_max=None):
    '''
    Easy set the limits of the axis.
    '''
    ax.set_xlim(left=x_min, right=x_max)
    ax.set_ylim(bottom=y_min, top=y_max)

def set_labels(ax, xlabel, ylabel):
    '''
    Easy set the labels of the axis.
    '''
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
