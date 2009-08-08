@set icc=D:\Program Files\Intel\Compiler\11.1\038\
@set sdk=D:\Program Files\Microsoft SDKs\Windows\v6.1
@set vs=D:\Program Files\Microsoft Visual Studio 9.0

@set lib=%icc%\lib\ia32;%icc%\tbb\ia32\vc9\lib;%vs%\VC\lib;%sdk%\Lib;.\
@set include=.\src;%icc%\include;%vs%\VC\include;%sdk%\Include;.\tbb\src;.\tbb\include;.\tbb\src;.\png
@set path=%icc%\bin\ia32;%vs%\Common7\IDE;%vs%\VC\bin;%sdk%\bin

@set librs=tbb-static.lib User32.Lib Kernel32.Lib Gdi32.Lib Winmm.lib png.lib
@set linkopt=ray.res /NODEFAULTLIB:tbb.lib /MACHINE:X86
@set defins=/DUSE_WINTHREAD /D__TBB_TASK_CPP_DIRECTLY_INCLUDED=1 /D_WIN32 /DUNICODE