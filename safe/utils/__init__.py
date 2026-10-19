from contextlib import contextmanager


@contextmanager
def attr_as(obj, field, value):
    """临时替换对象属性，退出上下文时恢复原值。"""
    old_value = getattr(obj, field, None)
    setattr(obj, field, value)
    try:
        yield
    finally:
        setattr(obj, field, old_value)