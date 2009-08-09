del theRayTracingGameRu.exe
mkdir TheRayTracingGameRu
copy distr-ru\*.* TheRayTracingGameRu
7z a -mx9 -sfx7z.sfx theRayTracingGameRu.exe TheRayTracingGameRu
del /Q TheRayTracingGameRu\*.* 
rmdir TheRayTracingGameRu
upx -9 theRayTracingGameRu.exe

del theRayTracingGameEn.exe
mkdir TheRayTracingGameEn
copy distr-en\*.* TheRayTracingGameEn
7z a -mx9 -sfx7z.sfx theRayTracingGameEn.exe TheRayTracingGameEn
del /Q TheRayTracingGameEn\*.* 
rmdir TheRayTracingGameEn
upx -9 theRayTracingGameEn.exe
