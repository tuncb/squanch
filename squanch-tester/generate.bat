@echo off

conan install . --install-folder build || goto :error
conan build . --build-folder build || goto :error

goto :success

:error
echo Failed with error #%errorlevel%.
exit /b %errorlevel%

:success
echo Success!
