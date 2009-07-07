@call env.bat
@rem /Wall /Wcheck /Wport /FAs /O3

set librs=tbb-static.lib User32.Lib Kernel32.Lib Gdi32.Lib Winmm.lib
set linkopt=ray.res /NODEFAULTLIB:tbb.lib /MACHINE:X86
set defins=/DUSE_WINTHREAD /D__TBB_TASK_CPP_DIRECTLY_INCLUDED=1 /D_WIN32 /DUNICODE
@rem /O3
icl /GX /fast /Qprof-gen %defins% /nologo /MT /Fo.\ /Fe.\distr\ray.exe .\src\main.cpp %librs% /link %linkopt% >build.log

@rem copy build.log con
