#Compass
## lib

The folder stores external libraries referenced by Compass which are not provided by Compass.

### standard references
Copy reference DLL and XML files to this directory.
You need to use the RhinoCommon version (not Grasshopper)
that comes with Rhino5, and is to find inside the
Rhinoceros 5.0\System folder

- GH_IO.dll
- GH_IO.xml
- Grasshopper.dll
- Grasshopper.xml
- RhinoCommon.dll
- RhinoCommon.xml


### other references
- Licenses for each external library should be named with the same filename and the `.license` extension.  (Not relevant to Grasshopper.dll, GH_IO.dll or RhinoCommons.dll)
e.g.:

    anexternallibrary.1.X.X.dll
    anexternallibrary.1.X.X.license
