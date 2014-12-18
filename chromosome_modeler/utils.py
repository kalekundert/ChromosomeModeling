from itertools import *

def pairwise(iterable):
    "s -> (s0,s1), (s1,s2), (s2, s3), ..."
    a, b = tee(iterable)
    next(b, None)
    return izip(a, b)

def check_all_kwargs_used(kwargs):
    kwargs = kwargs.copy()
    if kwargs and kwargs.pop('check_args', True):
        message = "unrecognized keyword argument{}: {}"
        plural = '' if len(kwargs) == 1 else 's'
        unknown_args = ', '.join("'{}'".format(k) for k in kwargs)
        raise TypeError(message.format(plural, unknown_args))



