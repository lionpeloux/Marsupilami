.. Marsupilami documentation master file, created by
   sphinx-quickstart on Sun Oct 18 13:44:27 2015.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

*******************************
API autodoc (WIP)
*******************************

Test handwritten doc (sphinxcontrib-dotnetdomain)
---------------------------------

This should be available in a near futur with autodoc capabilities with the `sphinx-autoapi <https://github.com/rtfd/sphinx-autoapi>`_ (WIP).

  .. dn:method:: Foo.method(FooType bar)

    :param bar: A Bar instance

    :type bar: :dn:cls:`Bar`

    :returns: Altered bar instance

  On :dn:cls:`Foo` instance, return :dn:cls:`Bar` instance

  .. dn:class:: ValidClass
    :protected:
    :static:

    .. dn:method:: MethodNoArgs()

    .. dn:method:: MethodArg(T1)

      :param T1: Version slug to use for node lookup

    .. dn:method:: MethodArgs(T1, T2)

      :param T1: desciption T1
      :param T2: desciption T2

    .. dn:method:: Foo.method(FooType bar)

      :param bar: A Bar instance
      :type bar: :dn:cls:`Bar`
      :returns: Altered bar instance

    .. dn:method:: MethodNested(List<int>, Dictionary<string,List<int>>)

    .. dn:property:: Foobar()

      :setter:
      :getter:

    .. dn:field:: Foobar()

      :adder:
      :remover:

  .. dn:namespace:: ValidNamespace

  .. dn:class:: ValidClass

  .. dn:struct:: ValidStructure

  .. dn:interface:: ValidInterface

  .. dn:property:: ValidProperty

  .. dn:field:: ValidField

  .. dn:event:: ValidEvent

  .. dn:operator:: ValidOperator

  .. dn:namespace:: ValidNamespace

    .. dn:class:: Foobar<T>

    .. dn:class:: Foobar<T,T>

    .. dn:class:: Foobar<TFoo,TBar>

    .. dn:class:: Foobar<T,<string>>

    .. dn:class:: Foobar<T,<T,<string>>>

    .. dn:property:: NestedProperty

    .. dn:field:: NestedField

    .. dn:event:: NestedEvent

    .. dn:operator:: NestedOperator

    .. dn:class:: NestedClass

      .. dn:property:: NestedClassProperty

      .. dn:field:: NestedClassField

      .. dn:event:: NestedClassEvent

      .. dn:operator:: NestedClassOperator


Dynamic Relaxation (breathe)
----------------------------

  .. doxygenclass:: Marsupilami::Kernel::DRRelax

Element (breathe)
-----------------

  .. doxygenclass:: Marsupilami::Kernel::DRElement

Utility (breathe)
-----------------

  .. doxygenclass:: Marsupilami::Kernel::Utility
