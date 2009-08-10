@call env.bat
@call generate-sources.bat

del pgopti.dpi.lock
del pgopti.dpi

icl /GX /fast /Qprof-use %defins% /nologo /MT /FAs /Fo.\ /Fe.\distr-ru\ray.pgo.exe .\src\main.cpp %librs% /link %linkopt% >build-pgo2.log