del theRayTracingGame.exe
mkdir TheRayTracingGame
copy distr\*.* TheRayTracingGame
7z a -sfx7z.sfx theRayTracingGame.exe TheRayTracingGame
del /Q TheRayTracingGame\*.* 
rmdir TheRayTracingGame
upx -9 theRayTracingGame.exe

