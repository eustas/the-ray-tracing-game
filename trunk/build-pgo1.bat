@call env.bat
@call generate-sources.bat

del *.dyn

icl /GX /fast /Qprof-gen %defins% /nologo /MT /Fo.\ /Fe.\distr-ru\ray.dyn.exe .\src\main.cpp %librs% /link %linkopt% >build-pgo1.log