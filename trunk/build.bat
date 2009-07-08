@call env.bat
@rem /Wall /Wcheck /Wport /FAs /O3

del  src\render-by-lines.class.Aa.cpp
copy src\render-by-lines.class.cpp src\render-by-lines.class.Aa.cpp
del  src\render-by-lines.class.Gp.cpp
copy src\render-by-lines.class.cpp src\render-by-lines.class.Gp.cpp
del  src\render-by-lines.class.AaGp.cpp
copy src\render-by-lines.class.cpp src\render-by-lines.class.AaGp.cpp

set librs=tbb-static.lib User32.Lib Kernel32.Lib Gdi32.Lib Winmm.lib png.lib
set linkopt=ray.res /NODEFAULTLIB:tbb.lib /MACHINE:X86
set defins=/DUSE_WINTHREAD /D__TBB_TASK_CPP_DIRECTLY_INCLUDED=1 /D_WIN32 /DUNICODE
@rem /O3 /P /Qipo
icl /GX /fast %defins% /nologo /MT /Fo.\ /Fe.\distr\ray.exe .\src\main.cpp %librs% /link %linkopt% >build.log

@rem copy build.log con
