# Microsoft Developer Studio Project File - Name="camb" - Package Owner=<4>
# Microsoft Developer Studio Generated Build File, Format Version 6.00
# ** DO NOT EDIT **

# TARGTYPE "Win32 (x86) Console Application" 0x0103

CFG=camb - Win32 Debug
!MESSAGE This is not a valid makefile. To build this project using NMAKE,
!MESSAGE use the Export Makefile command and run
!MESSAGE 
!MESSAGE NMAKE /f "camb.mak".
!MESSAGE 
!MESSAGE You can specify a configuration when running NMAKE
!MESSAGE by defining the macro CFG on the command line. For example:
!MESSAGE 
!MESSAGE NMAKE /f "camb.mak" CFG="camb - Win32 Debug"
!MESSAGE 
!MESSAGE Possible choices for configuration are:
!MESSAGE 
!MESSAGE "camb - Win32 Release" (based on "Win32 (x86) Console Application")
!MESSAGE "camb - Win32 Debug" (based on "Win32 (x86) Console Application")
!MESSAGE 

# Begin Project
# PROP AllowPerConfigDependencies 0
# PROP Scc_ProjName ""
# PROP Scc_LocalPath ""
CPP=cl.exe
F90=df.exe
RSC=rc.exe

!IF  "$(CFG)" == "camb - Win32 Release"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 0
# PROP BASE Output_Dir "Release"
# PROP BASE Intermediate_Dir "Release"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 0
# PROP Output_Dir "Release"
# PROP Intermediate_Dir "Release"
# PROP Ignore_Export_Lib 0
# PROP Target_Dir ""
# ADD BASE F90 /compile_only /nologo /warn:nofileopt
# ADD F90 /compile_only /fpp /nologo /warn:nofileopt
# ADD BASE CPP /nologo /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_CONSOLE" /D "_MBCS" /YX /FD /c
# ADD CPP /nologo /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_CONSOLE" /D "_MBCS" /YX /FD /c
# ADD BASE RSC /l 0x809 /d "NDEBUG"
# ADD RSC /l 0x809 /d "NDEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LINK32=link.exe
# ADD BASE LINK32 kernel32.lib /nologo /subsystem:console /machine:I386
# ADD LINK32 kernel32.lib /nologo /subsystem:console /machine:I386

!ELSEIF  "$(CFG)" == "camb - Win32 Debug"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 1
# PROP BASE Output_Dir "Debug"
# PROP BASE Intermediate_Dir "Debug"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 1
# PROP Output_Dir "Debug"
# PROP Intermediate_Dir "Debug"
# PROP Ignore_Export_Lib 0
# PROP Target_Dir ""
# ADD BASE F90 /check:bounds /compile_only /dbglibs /debug:full /nologo /traceback /warn:argument_checking /warn:nofileopt
# ADD F90 /check:bounds /compile_only /dbglibs /debug:full /fpp /nologo /traceback /warn:argument_checking /warn:nofileopt
# ADD BASE CPP /nologo /W3 /Gm /GX /ZI /Od /D "WIN32" /D "_DEBUG" /D "_CONSOLE" /D "_MBCS" /YX /FD /GZ /c
# ADD CPP /nologo /W3 /Gm /GX /ZI /Od /D "WIN32" /D "_DEBUG" /D "_CONSOLE" /D "_MBCS" /YX /FD /GZ /c
# ADD BASE RSC /l 0x809 /d "_DEBUG"
# ADD RSC /l 0x809 /d "_DEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LINK32=link.exe
# ADD BASE LINK32 kernel32.lib /nologo /subsystem:console /debug /machine:I386 /pdbtype:sept
# ADD LINK32 kernel32.lib /nologo /subsystem:console /incremental:no /debug /machine:I386 /out:"camb.exe" /pdbtype:sept

!ENDIF 

# Begin Target

# Name "camb - Win32 Release"
# Name "camb - Win32 Debug"
# Begin Group "Source Files"

# PROP Default_Filter "cpp;c;cxx;rc;def;r;odl;idl;hpj;bat;f90;for;f;fpp"
# Begin Source File

SOURCE=.\bessels.f90
NODEP_F90_BESSE=\
	".\Debug\lvalues.mod"\
	".\Debug\ModelParams.mod"\
	".\Debug\Precision.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\camb.f90
NODEP_F90_CAMB_=\
	".\Debug\CAMBmain.mod"\
	".\Debug\GaugeInterface.mod"\
	".\Debug\InitialPower.mod"\
	".\Debug\lensing.mod"\
	".\Debug\ModelData.mod"\
	".\Debug\ModelParams.mod"\
	".\Debug\Precision.mod"\
	".\Debug\Transfer.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\cmbmain.f90
NODEP_F90_CMBMA=\
	".\Debug\GaugeInterface.mod"\
	".\Debug\InitialPower.mod"\
	".\Debug\lvalues.mod"\
	".\Debug\MassiveNu.mod"\
	".\Debug\ModelData.mod"\
	".\Debug\ModelParams.mod"\
	".\Debug\NonLinear.mod"\
	".\Debug\Precision.mod"\
	".\Debug\SpherBessels.mod"\
	".\Debug\ThermoData.mod"\
	".\Debug\TimeSteps.mod"\
	".\Debug\Transfer.mod"\
	".\Debug\USpherBessels.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\constants.f90
# End Source File
# Begin Source File

SOURCE=.\equations.f90
NODEP_F90_EQUAT=\
	".\Debug\lvalues.mod"\
	".\Debug\MassiveNu.mod"\
	".\Debug\ModelData.mod"\
	".\Debug\ModelParams.mod"\
	".\Debug\Precision.mod"\
	".\Debug\ThermoData.mod"\
	".\Debug\Transfer.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\halofit.f90
NODEP_F90_HALOF=\
	".\Debug\ModelParams.mod"\
	".\Debug\Transfer.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\inidriver.F90
NODEP_F90_INIDR=\
	".\Debug\AmlUtils.mod"\
	".\Debug\CAMB.mod"\
	".\Debug\F90_UNIX.mod"\
	".\Debug\IniFile.mod"\
	".\Debug\LambdaGeneral.mod"\
	".\Debug\lensing.mod"\
	".\Debug\RECFAST.MOD"\
	
# End Source File
# Begin Source File

SOURCE=.\inifile.f90
# End Source File
# Begin Source File

SOURCE=.\lensing.f90
NODEP_F90_LENSI=\
	".\Debug\InitialPower.mod"\
	".\Debug\lvalues.mod"\
	".\Debug\ModelData.mod"\
	".\Debug\ModelParams.mod"\
	".\Debug\Precision.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\modules.f90
NODEP_F90_MODUL=\
	".\Debug\AmlUtils.mod"\
	".\Debug\IniFile.mod"\
	".\Debug\InitialPower.mod"\
	".\Debug\Precision.mod"\
	".\Debug\RECFAST.MOD"\
	
# End Source File
# Begin Source File

SOURCE=.\power_tilt.f90
NODEP_F90_POWER=\
	".\Debug\Precision.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\recfast.f90
NODEP_F90_RECFA=\
	".\Debug\Precision.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\reionization.f90
# End Source File
# Begin Source File

SOURCE=.\subroutines.f90
# End Source File
# Begin Source File

SOURCE=.\utils.F90
NODEP_F90_UTILS=\
	".\Debug\F90_UNIX.mod"\
	".\Debug\IFPORT.mod"\
	".\Debug\mpif.h"\
	".\xml_dll_use.mod"\
	".\XML_INCLUDE.F90"\
	".\xml_static_use.mod"\
	
# End Source File
# End Group
# Begin Group "Header Files"

# PROP Default_Filter "h;hpp;hxx;hm;inl;fi;fd"
# End Group
# Begin Group "Resource Files"

# PROP Default_Filter "ico;cur;bmp;dlg;rc2;rct;bin;rgs;gif;jpg;jpeg;jpe"
# End Group
# End Target
# End Project
