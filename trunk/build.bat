@call env.bat
@call generate-sources.bat

icl /GX %defins% /fast /nologo /MT /Fo.\ /Fe.\distr-ru\ray.exe .\src\main.cpp %librs% /link %linkopt% >build.log