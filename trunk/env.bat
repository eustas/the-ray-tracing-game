@set icc=D:\Program Files\Intel\Compiler\11.0\072\cpp
@set sdk=D:\Program Files\Microsoft SDKs\Windows\v6.0A
@set vs=D:\Program Files\Microsoft Visual Studio 9.0

@set lib=%icc%\lib\ia32;%icc%\tbb\ia32\vc9\lib;%vs%\VC\lib;%sdk%\Lib;.\
@set include=.\src;%icc%\include;%vs%\VC\include;%sdk%\Include;.\tbb\src;.\tbb\include;.\tbb\src
@set path=%icc%\bin\ia32;%vs%\Common7\IDE;%vs%\VC\bin;%sdk%\bin
