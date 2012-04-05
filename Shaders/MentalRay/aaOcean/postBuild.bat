SET XSI_INSTALL_DIR=C:\Users\amaan\Autodesk\Softimage_Winchester_RC1\Application

copy spdls\aaOceanDisplaceShader.spdl "%XSI_INSTALL_DIR%"\spdl
copy spdls\aaOceanFoamShader.spdl "%XSI_INSTALL_DIR%"\spdl
copy spdls\aaOceanImgVecDisplace.spdl "%XSI_INSTALL_DIR%"\spdl
REM copy spdls\aaOceanNormalsShader.spdl "%XSI_INSTALL_DIR%"\spdl

copy ..\Output\aaOceanMentalRay.dll "%XSI_INSTALL_DIR%"\bin\nt-x86-64

del ..\Output\aaOceanMentalRay.exp
del ..\Output\aaOceanMentalRay.lib
rd /S /Q x64