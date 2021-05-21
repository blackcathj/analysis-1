### Example Analysis Skeleton

This package is intended to serve as a minimal working example of an analysis package. You may use it to develop your own analysis, or, much more preferably, you should use it as a template to create your own named analysis package (e.g. DirectPhotonAnalysis or whatever is appropriate).

To do so, just search each file in this package for `MyNewAnalysis` and replace it with your desired package name. Take note of the different cases, e.g. the class name is MyNewAnalysis but there are some instances of `mynewanalysis` in the Makefile.

You may also use the already created perl script from Chris Pinkenburg to create your analysis package. Just run at the command line:

```
CreateSubsysRecoModule.pl <Module Name> --all
```

where `<Module Name>` is replaced by your desired name, e.g. `DirectPhotonAnalysis`.
