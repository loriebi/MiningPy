# This file was automatically generated by SWIG (http://www.swig.org).
# Version 2.0.10
#
# Do not make changes to this file unless you know what you are doing--modify
# the SWIG interface file instead.



from sys import version_info
if version_info >= (2,6,0):
    def swig_import_helper():
        from os.path import dirname
        import imp
        fp = None
        try:
            fp, pathname, description = imp.find_module('_averager', [dirname(__file__)])
        except ImportError:
            import _averager
            return _averager
        if fp is not None:
            try:
                _mod = imp.load_module('_averager', fp, pathname, description)
            finally:
                fp.close()
            return _mod
    _averager = swig_import_helper()
    del swig_import_helper
else:
    import _averager
del version_info
try:
    _swig_property = property
except NameError:
    pass # Python < 2.2 doesn't have 'property'.
def _swig_setattr_nondynamic(self,class_type,name,value,static=1):
    if (name == "thisown"): return self.this.own(value)
    if (name == "this"):
        if type(value).__name__ == 'SwigPyObject':
            self.__dict__[name] = value
            return
    method = class_type.__swig_setmethods__.get(name,None)
    if method: return method(self,value)
    if (not static):
        self.__dict__[name] = value
    else:
        raise AttributeError("You cannot add attributes to %s" % self)

def _swig_setattr(self,class_type,name,value):
    return _swig_setattr_nondynamic(self,class_type,name,value,0)

def _swig_getattr(self,class_type,name):
    if (name == "thisown"): return self.this.own()
    method = class_type.__swig_getmethods__.get(name,None)
    if method: return method(self)
    raise AttributeError(name)

def _swig_repr(self):
    try: strthis = "proxy of " + self.this.__repr__()
    except: strthis = ""
    return "<%s.%s; %s >" % (self.__class__.__module__, self.__class__.__name__, strthis,)

try:
    _object = object
    _newclass = 1
except AttributeError:
    class _object : pass
    _newclass = 0


class SwigPyIterator(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, SwigPyIterator, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, SwigPyIterator, name)
    def __init__(self, *args, **kwargs): raise AttributeError("No constructor defined - class is abstract")
    __repr__ = _swig_repr
    __swig_destroy__ = _averager.delete_SwigPyIterator
    __del__ = lambda self : None;
    def value(self): return _averager.SwigPyIterator_value(self)
    def incr(self, n=1): return _averager.SwigPyIterator_incr(self, n)
    def decr(self, n=1): return _averager.SwigPyIterator_decr(self, n)
    def distance(self, *args): return _averager.SwigPyIterator_distance(self, *args)
    def equal(self, *args): return _averager.SwigPyIterator_equal(self, *args)
    def copy(self): return _averager.SwigPyIterator_copy(self)
    def next(self): return _averager.SwigPyIterator_next(self)
    def __next__(self): return _averager.SwigPyIterator___next__(self)
    def previous(self): return _averager.SwigPyIterator_previous(self)
    def advance(self, *args): return _averager.SwigPyIterator_advance(self, *args)
    def __eq__(self, *args): return _averager.SwigPyIterator___eq__(self, *args)
    def __ne__(self, *args): return _averager.SwigPyIterator___ne__(self, *args)
    def __iadd__(self, *args): return _averager.SwigPyIterator___iadd__(self, *args)
    def __isub__(self, *args): return _averager.SwigPyIterator___isub__(self, *args)
    def __add__(self, *args): return _averager.SwigPyIterator___add__(self, *args)
    def __sub__(self, *args): return _averager.SwigPyIterator___sub__(self, *args)
    def __iter__(self): return self
SwigPyIterator_swigregister = _averager.SwigPyIterator_swigregister
SwigPyIterator_swigregister(SwigPyIterator)

class vecFl(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, vecFl, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, vecFl, name)
    __repr__ = _swig_repr
    def iterator(self): return _averager.vecFl_iterator(self)
    def __iter__(self): return self.iterator()
    def __nonzero__(self): return _averager.vecFl___nonzero__(self)
    def __bool__(self): return _averager.vecFl___bool__(self)
    def __len__(self): return _averager.vecFl___len__(self)
    def pop(self): return _averager.vecFl_pop(self)
    def __getslice__(self, *args): return _averager.vecFl___getslice__(self, *args)
    def __setslice__(self, *args): return _averager.vecFl___setslice__(self, *args)
    def __delslice__(self, *args): return _averager.vecFl___delslice__(self, *args)
    def __delitem__(self, *args): return _averager.vecFl___delitem__(self, *args)
    def __getitem__(self, *args): return _averager.vecFl___getitem__(self, *args)
    def __setitem__(self, *args): return _averager.vecFl___setitem__(self, *args)
    def append(self, *args): return _averager.vecFl_append(self, *args)
    def empty(self): return _averager.vecFl_empty(self)
    def size(self): return _averager.vecFl_size(self)
    def clear(self): return _averager.vecFl_clear(self)
    def swap(self, *args): return _averager.vecFl_swap(self, *args)
    def get_allocator(self): return _averager.vecFl_get_allocator(self)
    def begin(self): return _averager.vecFl_begin(self)
    def end(self): return _averager.vecFl_end(self)
    def rbegin(self): return _averager.vecFl_rbegin(self)
    def rend(self): return _averager.vecFl_rend(self)
    def pop_back(self): return _averager.vecFl_pop_back(self)
    def erase(self, *args): return _averager.vecFl_erase(self, *args)
    def __init__(self, *args): 
        this = _averager.new_vecFl(*args)
        try: self.this.append(this)
        except: self.this = this
    def push_back(self, *args): return _averager.vecFl_push_back(self, *args)
    def front(self): return _averager.vecFl_front(self)
    def back(self): return _averager.vecFl_back(self)
    def assign(self, *args): return _averager.vecFl_assign(self, *args)
    def resize(self, *args): return _averager.vecFl_resize(self, *args)
    def insert(self, *args): return _averager.vecFl_insert(self, *args)
    def reserve(self, *args): return _averager.vecFl_reserve(self, *args)
    def capacity(self): return _averager.vecFl_capacity(self)
    __swig_destroy__ = _averager.delete_vecFl
    __del__ = lambda self : None;
vecFl_swigregister = _averager.vecFl_swigregister
vecFl_swigregister(vecFl)

class vecDb(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, vecDb, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, vecDb, name)
    __repr__ = _swig_repr
    def iterator(self): return _averager.vecDb_iterator(self)
    def __iter__(self): return self.iterator()
    def __nonzero__(self): return _averager.vecDb___nonzero__(self)
    def __bool__(self): return _averager.vecDb___bool__(self)
    def __len__(self): return _averager.vecDb___len__(self)
    def pop(self): return _averager.vecDb_pop(self)
    def __getslice__(self, *args): return _averager.vecDb___getslice__(self, *args)
    def __setslice__(self, *args): return _averager.vecDb___setslice__(self, *args)
    def __delslice__(self, *args): return _averager.vecDb___delslice__(self, *args)
    def __delitem__(self, *args): return _averager.vecDb___delitem__(self, *args)
    def __getitem__(self, *args): return _averager.vecDb___getitem__(self, *args)
    def __setitem__(self, *args): return _averager.vecDb___setitem__(self, *args)
    def append(self, *args): return _averager.vecDb_append(self, *args)
    def empty(self): return _averager.vecDb_empty(self)
    def size(self): return _averager.vecDb_size(self)
    def clear(self): return _averager.vecDb_clear(self)
    def swap(self, *args): return _averager.vecDb_swap(self, *args)
    def get_allocator(self): return _averager.vecDb_get_allocator(self)
    def begin(self): return _averager.vecDb_begin(self)
    def end(self): return _averager.vecDb_end(self)
    def rbegin(self): return _averager.vecDb_rbegin(self)
    def rend(self): return _averager.vecDb_rend(self)
    def pop_back(self): return _averager.vecDb_pop_back(self)
    def erase(self, *args): return _averager.vecDb_erase(self, *args)
    def __init__(self, *args): 
        this = _averager.new_vecDb(*args)
        try: self.this.append(this)
        except: self.this = this
    def push_back(self, *args): return _averager.vecDb_push_back(self, *args)
    def front(self): return _averager.vecDb_front(self)
    def back(self): return _averager.vecDb_back(self)
    def assign(self, *args): return _averager.vecDb_assign(self, *args)
    def resize(self, *args): return _averager.vecDb_resize(self, *args)
    def insert(self, *args): return _averager.vecDb_insert(self, *args)
    def reserve(self, *args): return _averager.vecDb_reserve(self, *args)
    def capacity(self): return _averager.vecDb_capacity(self)
    __swig_destroy__ = _averager.delete_vecDb
    __del__ = lambda self : None;
vecDb_swigregister = _averager.vecDb_swigregister
vecDb_swigregister(vecDb)

class vecUint(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, vecUint, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, vecUint, name)
    __repr__ = _swig_repr
    def iterator(self): return _averager.vecUint_iterator(self)
    def __iter__(self): return self.iterator()
    def __nonzero__(self): return _averager.vecUint___nonzero__(self)
    def __bool__(self): return _averager.vecUint___bool__(self)
    def __len__(self): return _averager.vecUint___len__(self)
    def pop(self): return _averager.vecUint_pop(self)
    def __getslice__(self, *args): return _averager.vecUint___getslice__(self, *args)
    def __setslice__(self, *args): return _averager.vecUint___setslice__(self, *args)
    def __delslice__(self, *args): return _averager.vecUint___delslice__(self, *args)
    def __delitem__(self, *args): return _averager.vecUint___delitem__(self, *args)
    def __getitem__(self, *args): return _averager.vecUint___getitem__(self, *args)
    def __setitem__(self, *args): return _averager.vecUint___setitem__(self, *args)
    def append(self, *args): return _averager.vecUint_append(self, *args)
    def empty(self): return _averager.vecUint_empty(self)
    def size(self): return _averager.vecUint_size(self)
    def clear(self): return _averager.vecUint_clear(self)
    def swap(self, *args): return _averager.vecUint_swap(self, *args)
    def get_allocator(self): return _averager.vecUint_get_allocator(self)
    def begin(self): return _averager.vecUint_begin(self)
    def end(self): return _averager.vecUint_end(self)
    def rbegin(self): return _averager.vecUint_rbegin(self)
    def rend(self): return _averager.vecUint_rend(self)
    def pop_back(self): return _averager.vecUint_pop_back(self)
    def erase(self, *args): return _averager.vecUint_erase(self, *args)
    def __init__(self, *args): 
        this = _averager.new_vecUint(*args)
        try: self.this.append(this)
        except: self.this = this
    def push_back(self, *args): return _averager.vecUint_push_back(self, *args)
    def front(self): return _averager.vecUint_front(self)
    def back(self): return _averager.vecUint_back(self)
    def assign(self, *args): return _averager.vecUint_assign(self, *args)
    def resize(self, *args): return _averager.vecUint_resize(self, *args)
    def insert(self, *args): return _averager.vecUint_insert(self, *args)
    def reserve(self, *args): return _averager.vecUint_reserve(self, *args)
    def capacity(self): return _averager.vecUint_capacity(self)
    __swig_destroy__ = _averager.delete_vecUint
    __del__ = lambda self : None;
vecUint_swigregister = _averager.vecUint_swigregister
vecUint_swigregister(vecUint)

class ParamsThread(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, ParamsThread, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, ParamsThread, name)
    __repr__ = _swig_repr
    __swig_setmethods__["index"] = _averager.ParamsThread_index_set
    __swig_getmethods__["index"] = _averager.ParamsThread_index_get
    if _newclass:index = _swig_property(_averager.ParamsThread_index_get, _averager.ParamsThread_index_set)
    __swig_setmethods__["startY"] = _averager.ParamsThread_startY_set
    __swig_getmethods__["startY"] = _averager.ParamsThread_startY_get
    if _newclass:startY = _swig_property(_averager.ParamsThread_startY_get, _averager.ParamsThread_startY_set)
    __swig_setmethods__["endY"] = _averager.ParamsThread_endY_set
    __swig_getmethods__["endY"] = _averager.ParamsThread_endY_get
    if _newclass:endY = _swig_property(_averager.ParamsThread_endY_get, _averager.ParamsThread_endY_set)
    def __init__(self): 
        this = _averager.new_ParamsThread()
        try: self.this.append(this)
        except: self.this = this
    __swig_destroy__ = _averager.delete_ParamsThread
    __del__ = lambda self : None;
ParamsThread_swigregister = _averager.ParamsThread_swigregister
ParamsThread_swigregister(ParamsThread)

class Averager(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, Averager, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, Averager, name)
    __repr__ = _swig_repr
    def __init__(self, *args): 
        this = _averager.new_Averager(*args)
        try: self.this.append(this)
        except: self.this = this
    def setVariables(self, *args): return _averager.Averager_setVariables(self, *args)
    def readImage(self, *args): return _averager.Averager_readImage(self, *args)
    def getDataDimentions(self): return _averager.Averager_getDataDimentions(self)
    def runThreads(self, *args): return _averager.Averager_runThreads(self, *args)
    def countAverage(self): return _averager.Averager_countAverage(self)
    __swig_destroy__ = _averager.delete_Averager
    __del__ = lambda self : None;
Averager_swigregister = _averager.Averager_swigregister
Averager_swigregister(Averager)

# This file is compatible with both classic and new-style classes.


