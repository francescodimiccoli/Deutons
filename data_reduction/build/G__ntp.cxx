// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME G__ntp

/*******************************************************************/
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#define G__DICTIONARY
#include "RConfig.h"
#include "TClass.h"
#include "TDictAttributeMap.h"
#include "TInterpreter.h"
#include "TROOT.h"
#include "TBuffer.h"
#include "TMemberInspector.h"
#include "TInterpreter.h"
#include "TVirtualMutex.h"
#include "TError.h"

#ifndef G__ROOT
#define G__ROOT
#endif

#include "RtypesImp.h"
#include "TIsAProxy.h"
#include "TFileMergeInfo.h"
#include <algorithm>
#include "TCollectionProxyInfo.h"
/*******************************************************************/

#include "TDataMember.h"

// Since CINT ignores the std namespace, we need to do so in this file.
namespace std {} using namespace std;

// Header files passed as explicit arguments
#include "/afs/cern.ch/user/f/fdimicco/Work/Deutons/data_reduction/include/Ntp.h"
#include "/afs/cern.ch/user/f/fdimicco/Work/Deutons/data_reduction/include/VPSArchive.h"
#include "/afs/cern.ch/user/f/fdimicco/Work/Deutons/data_reduction/include/RichOccupancy.h"
#include "/afs/cern.ch/user/f/fdimicco/Work/Deutons/data_reduction/include/Classifier.h"
#include "/afs/cern.ch/user/f/fdimicco/Work/Deutons/data_reduction/include/RichBDT.h"
#include "/afs/cern.ch/user/f/fdimicco/Work/Deutons/data_reduction/include/MSpline.h"
#include "/afs/cern.ch/user/f/fdimicco/Work/Deutons/data_reduction/include/IsoPoly.h"
#include "/afs/cern.ch/user/f/fdimicco/Work/Deutons/data_reduction/include/PDF1.h"
#include "/afs/cern.ch/user/f/fdimicco/Work/Deutons/data_reduction/include/PDF2.h"
#include "/afs/cern.ch/user/f/fdimicco/Work/Deutons/data_reduction/include/PDF2DB.h"
#include "/afs/cern.ch/user/f/fdimicco/Work/Deutons/data_reduction/include/LikelihoodVar.h"
#include "/afs/cern.ch/user/f/fdimicco/Work/Deutons/data_reduction/include/Likelihood.h"

// Header files passed via #pragma extra_include

namespace ROOT {
   static void *new_RTIInfo(void *p = 0);
   static void *newArray_RTIInfo(Long_t size, void *p);
   static void delete_RTIInfo(void *p);
   static void deleteArray_RTIInfo(void *p);
   static void destruct_RTIInfo(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::RTIInfo*)
   {
      ::RTIInfo *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::RTIInfo >(0);
      static ::ROOT::TGenericClassInfo 
         instance("RTIInfo", ::RTIInfo::Class_Version(), "Ntp.h", 22,
                  typeid(::RTIInfo), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::RTIInfo::Dictionary, isa_proxy, 4,
                  sizeof(::RTIInfo) );
      instance.SetNew(&new_RTIInfo);
      instance.SetNewArray(&newArray_RTIInfo);
      instance.SetDelete(&delete_RTIInfo);
      instance.SetDeleteArray(&deleteArray_RTIInfo);
      instance.SetDestructor(&destruct_RTIInfo);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::RTIInfo*)
   {
      return GenerateInitInstanceLocal((::RTIInfo*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::RTIInfo*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_FileInfo(void *p = 0);
   static void *newArray_FileInfo(Long_t size, void *p);
   static void delete_FileInfo(void *p);
   static void deleteArray_FileInfo(void *p);
   static void destruct_FileInfo(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::FileInfo*)
   {
      ::FileInfo *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::FileInfo >(0);
      static ::ROOT::TGenericClassInfo 
         instance("FileInfo", ::FileInfo::Class_Version(), "Ntp.h", 77,
                  typeid(::FileInfo), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::FileInfo::Dictionary, isa_proxy, 4,
                  sizeof(::FileInfo) );
      instance.SetNew(&new_FileInfo);
      instance.SetNewArray(&newArray_FileInfo);
      instance.SetDelete(&delete_FileInfo);
      instance.SetDeleteArray(&deleteArray_FileInfo);
      instance.SetDestructor(&destruct_FileInfo);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::FileInfo*)
   {
      return GenerateInitInstanceLocal((::FileInfo*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::FileInfo*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_FileMCInfo(void *p = 0);
   static void *newArray_FileMCInfo(Long_t size, void *p);
   static void delete_FileMCInfo(void *p);
   static void deleteArray_FileMCInfo(void *p);
   static void destruct_FileMCInfo(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::FileMCInfo*)
   {
      ::FileMCInfo *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::FileMCInfo >(0);
      static ::ROOT::TGenericClassInfo 
         instance("FileMCInfo", ::FileMCInfo::Class_Version(), "Ntp.h", 95,
                  typeid(::FileMCInfo), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::FileMCInfo::Dictionary, isa_proxy, 4,
                  sizeof(::FileMCInfo) );
      instance.SetNew(&new_FileMCInfo);
      instance.SetNewArray(&newArray_FileMCInfo);
      instance.SetDelete(&delete_FileMCInfo);
      instance.SetDeleteArray(&deleteArray_FileMCInfo);
      instance.SetDestructor(&destruct_FileMCInfo);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::FileMCInfo*)
   {
      return GenerateInitInstanceLocal((::FileMCInfo*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::FileMCInfo*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_ProcInfo(void *p = 0);
   static void *newArray_ProcInfo(Long_t size, void *p);
   static void delete_ProcInfo(void *p);
   static void deleteArray_ProcInfo(void *p);
   static void destruct_ProcInfo(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::ProcInfo*)
   {
      ::ProcInfo *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::ProcInfo >(0);
      static ::ROOT::TGenericClassInfo 
         instance("ProcInfo", ::ProcInfo::Class_Version(), "Ntp.h", 127,
                  typeid(::ProcInfo), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::ProcInfo::Dictionary, isa_proxy, 4,
                  sizeof(::ProcInfo) );
      instance.SetNew(&new_ProcInfo);
      instance.SetNewArray(&newArray_ProcInfo);
      instance.SetDelete(&delete_ProcInfo);
      instance.SetDeleteArray(&deleteArray_ProcInfo);
      instance.SetDestructor(&destruct_ProcInfo);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::ProcInfo*)
   {
      return GenerateInitInstanceLocal((::ProcInfo*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::ProcInfo*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_NtpSHeader(void *p = 0);
   static void *newArray_NtpSHeader(Long_t size, void *p);
   static void delete_NtpSHeader(void *p);
   static void deleteArray_NtpSHeader(void *p);
   static void destruct_NtpSHeader(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::NtpSHeader*)
   {
      ::NtpSHeader *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::NtpSHeader >(0);
      static ::ROOT::TGenericClassInfo 
         instance("NtpSHeader", ::NtpSHeader::Class_Version(), "Ntp.h", 148,
                  typeid(::NtpSHeader), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::NtpSHeader::Dictionary, isa_proxy, 4,
                  sizeof(::NtpSHeader) );
      instance.SetNew(&new_NtpSHeader);
      instance.SetNewArray(&newArray_NtpSHeader);
      instance.SetDelete(&delete_NtpSHeader);
      instance.SetDeleteArray(&deleteArray_NtpSHeader);
      instance.SetDestructor(&destruct_NtpSHeader);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::NtpSHeader*)
   {
      return GenerateInitInstanceLocal((::NtpSHeader*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::NtpSHeader*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_NtpHeader(void *p = 0);
   static void *newArray_NtpHeader(Long_t size, void *p);
   static void delete_NtpHeader(void *p);
   static void deleteArray_NtpHeader(void *p);
   static void destruct_NtpHeader(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::NtpHeader*)
   {
      ::NtpHeader *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::NtpHeader >(0);
      static ::ROOT::TGenericClassInfo 
         instance("NtpHeader", ::NtpHeader::Class_Version(), "Ntp.h", 165,
                  typeid(::NtpHeader), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::NtpHeader::Dictionary, isa_proxy, 4,
                  sizeof(::NtpHeader) );
      instance.SetNew(&new_NtpHeader);
      instance.SetNewArray(&newArray_NtpHeader);
      instance.SetDelete(&delete_NtpHeader);
      instance.SetDeleteArray(&deleteArray_NtpHeader);
      instance.SetDestructor(&destruct_NtpHeader);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::NtpHeader*)
   {
      return GenerateInitInstanceLocal((::NtpHeader*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::NtpHeader*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_NtpMCHeader(void *p = 0);
   static void *newArray_NtpMCHeader(Long_t size, void *p);
   static void delete_NtpMCHeader(void *p);
   static void deleteArray_NtpMCHeader(void *p);
   static void destruct_NtpMCHeader(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::NtpMCHeader*)
   {
      ::NtpMCHeader *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::NtpMCHeader >(0);
      static ::ROOT::TGenericClassInfo 
         instance("NtpMCHeader", ::NtpMCHeader::Class_Version(), "Ntp.h", 242,
                  typeid(::NtpMCHeader), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::NtpMCHeader::Dictionary, isa_proxy, 4,
                  sizeof(::NtpMCHeader) );
      instance.SetNew(&new_NtpMCHeader);
      instance.SetNewArray(&newArray_NtpMCHeader);
      instance.SetDelete(&delete_NtpMCHeader);
      instance.SetDeleteArray(&deleteArray_NtpMCHeader);
      instance.SetDestructor(&destruct_NtpMCHeader);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::NtpMCHeader*)
   {
      return GenerateInitInstanceLocal((::NtpMCHeader*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::NtpMCHeader*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_NtpTrd(void *p = 0);
   static void *newArray_NtpTrd(Long_t size, void *p);
   static void delete_NtpTrd(void *p);
   static void deleteArray_NtpTrd(void *p);
   static void destruct_NtpTrd(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::NtpTrd*)
   {
      ::NtpTrd *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::NtpTrd >(0);
      static ::ROOT::TGenericClassInfo 
         instance("NtpTrd", ::NtpTrd::Class_Version(), "Ntp.h", 272,
                  typeid(::NtpTrd), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::NtpTrd::Dictionary, isa_proxy, 4,
                  sizeof(::NtpTrd) );
      instance.SetNew(&new_NtpTrd);
      instance.SetNewArray(&newArray_NtpTrd);
      instance.SetDelete(&delete_NtpTrd);
      instance.SetDeleteArray(&deleteArray_NtpTrd);
      instance.SetDestructor(&destruct_NtpTrd);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::NtpTrd*)
   {
      return GenerateInitInstanceLocal((::NtpTrd*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::NtpTrd*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_NtpTof(void *p = 0);
   static void *newArray_NtpTof(Long_t size, void *p);
   static void delete_NtpTof(void *p);
   static void deleteArray_NtpTof(void *p);
   static void destruct_NtpTof(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::NtpTof*)
   {
      ::NtpTof *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::NtpTof >(0);
      static ::ROOT::TGenericClassInfo 
         instance("NtpTof", ::NtpTof::Class_Version(), "Ntp.h", 334,
                  typeid(::NtpTof), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::NtpTof::Dictionary, isa_proxy, 4,
                  sizeof(::NtpTof) );
      instance.SetNew(&new_NtpTof);
      instance.SetNewArray(&newArray_NtpTof);
      instance.SetDelete(&delete_NtpTof);
      instance.SetDeleteArray(&deleteArray_NtpTof);
      instance.SetDestructor(&destruct_NtpTof);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::NtpTof*)
   {
      return GenerateInitInstanceLocal((::NtpTof*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::NtpTof*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_NtpTracker(void *p = 0);
   static void *newArray_NtpTracker(Long_t size, void *p);
   static void delete_NtpTracker(void *p);
   static void deleteArray_NtpTracker(void *p);
   static void destruct_NtpTracker(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::NtpTracker*)
   {
      ::NtpTracker *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::NtpTracker >(0);
      static ::ROOT::TGenericClassInfo 
         instance("NtpTracker", ::NtpTracker::Class_Version(), "Ntp.h", 404,
                  typeid(::NtpTracker), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::NtpTracker::Dictionary, isa_proxy, 4,
                  sizeof(::NtpTracker) );
      instance.SetNew(&new_NtpTracker);
      instance.SetNewArray(&newArray_NtpTracker);
      instance.SetDelete(&delete_NtpTracker);
      instance.SetDeleteArray(&deleteArray_NtpTracker);
      instance.SetDestructor(&destruct_NtpTracker);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::NtpTracker*)
   {
      return GenerateInitInstanceLocal((::NtpTracker*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::NtpTracker*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_NtpRich(void *p = 0);
   static void *newArray_NtpRich(Long_t size, void *p);
   static void delete_NtpRich(void *p);
   static void deleteArray_NtpRich(void *p);
   static void destruct_NtpRich(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::NtpRich*)
   {
      ::NtpRich *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::NtpRich >(0);
      static ::ROOT::TGenericClassInfo 
         instance("NtpRich", ::NtpRich::Class_Version(), "Ntp.h", 488,
                  typeid(::NtpRich), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::NtpRich::Dictionary, isa_proxy, 4,
                  sizeof(::NtpRich) );
      instance.SetNew(&new_NtpRich);
      instance.SetNewArray(&newArray_NtpRich);
      instance.SetDelete(&delete_NtpRich);
      instance.SetDeleteArray(&deleteArray_NtpRich);
      instance.SetDestructor(&destruct_NtpRich);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::NtpRich*)
   {
      return GenerateInitInstanceLocal((::NtpRich*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::NtpRich*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_NtpEcal(void *p = 0);
   static void *newArray_NtpEcal(Long_t size, void *p);
   static void delete_NtpEcal(void *p);
   static void deleteArray_NtpEcal(void *p);
   static void destruct_NtpEcal(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::NtpEcal*)
   {
      ::NtpEcal *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::NtpEcal >(0);
      static ::ROOT::TGenericClassInfo 
         instance("NtpEcal", ::NtpEcal::Class_Version(), "Ntp.h", 600,
                  typeid(::NtpEcal), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::NtpEcal::Dictionary, isa_proxy, 4,
                  sizeof(::NtpEcal) );
      instance.SetNew(&new_NtpEcal);
      instance.SetNewArray(&newArray_NtpEcal);
      instance.SetDelete(&delete_NtpEcal);
      instance.SetDeleteArray(&deleteArray_NtpEcal);
      instance.SetDestructor(&destruct_NtpEcal);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::NtpEcal*)
   {
      return GenerateInitInstanceLocal((::NtpEcal*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::NtpEcal*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_NtpAnti(void *p = 0);
   static void *newArray_NtpAnti(Long_t size, void *p);
   static void delete_NtpAnti(void *p);
   static void deleteArray_NtpAnti(void *p);
   static void destruct_NtpAnti(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::NtpAnti*)
   {
      ::NtpAnti *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::NtpAnti >(0);
      static ::ROOT::TGenericClassInfo 
         instance("NtpAnti", ::NtpAnti::Class_Version(), "Ntp.h", 631,
                  typeid(::NtpAnti), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::NtpAnti::Dictionary, isa_proxy, 4,
                  sizeof(::NtpAnti) );
      instance.SetNew(&new_NtpAnti);
      instance.SetNewArray(&newArray_NtpAnti);
      instance.SetDelete(&delete_NtpAnti);
      instance.SetDeleteArray(&deleteArray_NtpAnti);
      instance.SetDestructor(&destruct_NtpAnti);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::NtpAnti*)
   {
      return GenerateInitInstanceLocal((::NtpAnti*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::NtpAnti*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_NtpStandAlone(void *p = 0);
   static void *newArray_NtpStandAlone(Long_t size, void *p);
   static void delete_NtpStandAlone(void *p);
   static void deleteArray_NtpStandAlone(void *p);
   static void destruct_NtpStandAlone(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::NtpStandAlone*)
   {
      ::NtpStandAlone *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::NtpStandAlone >(0);
      static ::ROOT::TGenericClassInfo 
         instance("NtpStandAlone", ::NtpStandAlone::Class_Version(), "Ntp.h", 659,
                  typeid(::NtpStandAlone), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::NtpStandAlone::Dictionary, isa_proxy, 4,
                  sizeof(::NtpStandAlone) );
      instance.SetNew(&new_NtpStandAlone);
      instance.SetNewArray(&newArray_NtpStandAlone);
      instance.SetDelete(&delete_NtpStandAlone);
      instance.SetDeleteArray(&deleteArray_NtpStandAlone);
      instance.SetDestructor(&destruct_NtpStandAlone);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::NtpStandAlone*)
   {
      return GenerateInitInstanceLocal((::NtpStandAlone*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::NtpStandAlone*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_Event(void *p = 0);
   static void *newArray_Event(Long_t size, void *p);
   static void delete_Event(void *p);
   static void deleteArray_Event(void *p);
   static void destruct_Event(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::Event*)
   {
      ::Event *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::Event >(0);
      static ::ROOT::TGenericClassInfo 
         instance("Event", ::Event::Class_Version(), "Ntp.h", 709,
                  typeid(::Event), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::Event::Dictionary, isa_proxy, 4,
                  sizeof(::Event) );
      instance.SetNew(&new_Event);
      instance.SetNewArray(&newArray_Event);
      instance.SetDelete(&delete_Event);
      instance.SetDeleteArray(&deleteArray_Event);
      instance.SetDestructor(&destruct_Event);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::Event*)
   {
      return GenerateInitInstanceLocal((::Event*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::Event*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_NtpCompact(void *p = 0);
   static void *newArray_NtpCompact(Long_t size, void *p);
   static void delete_NtpCompact(void *p);
   static void deleteArray_NtpCompact(void *p);
   static void destruct_NtpCompact(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::NtpCompact*)
   {
      ::NtpCompact *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::NtpCompact >(0);
      static ::ROOT::TGenericClassInfo 
         instance("NtpCompact", ::NtpCompact::Class_Version(), "Ntp.h", 756,
                  typeid(::NtpCompact), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::NtpCompact::Dictionary, isa_proxy, 4,
                  sizeof(::NtpCompact) );
      instance.SetNew(&new_NtpCompact);
      instance.SetNewArray(&newArray_NtpCompact);
      instance.SetDelete(&delete_NtpCompact);
      instance.SetDeleteArray(&deleteArray_NtpCompact);
      instance.SetDestructor(&destruct_NtpCompact);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::NtpCompact*)
   {
      return GenerateInitInstanceLocal((::NtpCompact*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::NtpCompact*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_VPSCategory(void *p = 0);
   static void *newArray_VPSCategory(Long_t size, void *p);
   static void delete_VPSCategory(void *p);
   static void deleteArray_VPSCategory(void *p);
   static void destruct_VPSCategory(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::VPSCategory*)
   {
      ::VPSCategory *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::VPSCategory >(0);
      static ::ROOT::TGenericClassInfo 
         instance("VPSCategory", ::VPSCategory::Class_Version(), "VPSArchive.h", 17,
                  typeid(::VPSCategory), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::VPSCategory::Dictionary, isa_proxy, 4,
                  sizeof(::VPSCategory) );
      instance.SetNew(&new_VPSCategory);
      instance.SetNewArray(&newArray_VPSCategory);
      instance.SetDelete(&delete_VPSCategory);
      instance.SetDeleteArray(&deleteArray_VPSCategory);
      instance.SetDestructor(&destruct_VPSCategory);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::VPSCategory*)
   {
      return GenerateInitInstanceLocal((::VPSCategory*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::VPSCategory*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_VPSArchive(void *p = 0);
   static void *newArray_VPSArchive(Long_t size, void *p);
   static void delete_VPSArchive(void *p);
   static void deleteArray_VPSArchive(void *p);
   static void destruct_VPSArchive(void *p);
   static Long64_t merge_VPSArchive(void *obj, TCollection *coll,TFileMergeInfo *info);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::VPSArchive*)
   {
      ::VPSArchive *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::VPSArchive >(0);
      static ::ROOT::TGenericClassInfo 
         instance("VPSArchive", ::VPSArchive::Class_Version(), "VPSArchive.h", 82,
                  typeid(::VPSArchive), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::VPSArchive::Dictionary, isa_proxy, 4,
                  sizeof(::VPSArchive) );
      instance.SetNew(&new_VPSArchive);
      instance.SetNewArray(&newArray_VPSArchive);
      instance.SetDelete(&delete_VPSArchive);
      instance.SetDeleteArray(&deleteArray_VPSArchive);
      instance.SetDestructor(&destruct_VPSArchive);
      instance.SetMerge(&merge_VPSArchive);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::VPSArchive*)
   {
      return GenerateInitInstanceLocal((::VPSArchive*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::VPSArchive*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_RichOccupancy(void *p = 0);
   static void *newArray_RichOccupancy(Long_t size, void *p);
   static void delete_RichOccupancy(void *p);
   static void deleteArray_RichOccupancy(void *p);
   static void destruct_RichOccupancy(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::RichOccupancy*)
   {
      ::RichOccupancy *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::RichOccupancy >(0);
      static ::ROOT::TGenericClassInfo 
         instance("RichOccupancy", ::RichOccupancy::Class_Version(), "RichOccupancy.h", 17,
                  typeid(::RichOccupancy), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::RichOccupancy::Dictionary, isa_proxy, 4,
                  sizeof(::RichOccupancy) );
      instance.SetNew(&new_RichOccupancy);
      instance.SetNewArray(&newArray_RichOccupancy);
      instance.SetDelete(&delete_RichOccupancy);
      instance.SetDeleteArray(&deleteArray_RichOccupancy);
      instance.SetDestructor(&destruct_RichOccupancy);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::RichOccupancy*)
   {
      return GenerateInitInstanceLocal((::RichOccupancy*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::RichOccupancy*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_IsoPoly(void *p = 0);
   static void *newArray_IsoPoly(Long_t size, void *p);
   static void delete_IsoPoly(void *p);
   static void deleteArray_IsoPoly(void *p);
   static void destruct_IsoPoly(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::IsoPoly*)
   {
      ::IsoPoly *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::IsoPoly >(0);
      static ::ROOT::TGenericClassInfo 
         instance("IsoPoly", ::IsoPoly::Class_Version(), "IsoPoly.h", 10,
                  typeid(::IsoPoly), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::IsoPoly::Dictionary, isa_proxy, 4,
                  sizeof(::IsoPoly) );
      instance.SetNew(&new_IsoPoly);
      instance.SetNewArray(&newArray_IsoPoly);
      instance.SetDelete(&delete_IsoPoly);
      instance.SetDeleteArray(&deleteArray_IsoPoly);
      instance.SetDestructor(&destruct_IsoPoly);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::IsoPoly*)
   {
      return GenerateInitInstanceLocal((::IsoPoly*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::IsoPoly*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_MSpline(void *p = 0);
   static void *newArray_MSpline(Long_t size, void *p);
   static void delete_MSpline(void *p);
   static void deleteArray_MSpline(void *p);
   static void destruct_MSpline(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::MSpline*)
   {
      ::MSpline *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::MSpline >(0);
      static ::ROOT::TGenericClassInfo 
         instance("MSpline", ::MSpline::Class_Version(), "MSpline.h", 15,
                  typeid(::MSpline), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::MSpline::Dictionary, isa_proxy, 4,
                  sizeof(::MSpline) );
      instance.SetNew(&new_MSpline);
      instance.SetNewArray(&newArray_MSpline);
      instance.SetDelete(&delete_MSpline);
      instance.SetDeleteArray(&deleteArray_MSpline);
      instance.SetDestructor(&destruct_MSpline);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::MSpline*)
   {
      return GenerateInitInstanceLocal((::MSpline*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::MSpline*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_PDF1(void *p = 0);
   static void *newArray_PDF1(Long_t size, void *p);
   static void delete_PDF1(void *p);
   static void deleteArray_PDF1(void *p);
   static void destruct_PDF1(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::PDF1*)
   {
      ::PDF1 *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::PDF1 >(0);
      static ::ROOT::TGenericClassInfo 
         instance("PDF1", ::PDF1::Class_Version(), "PDF1.h", 20,
                  typeid(::PDF1), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::PDF1::Dictionary, isa_proxy, 4,
                  sizeof(::PDF1) );
      instance.SetNew(&new_PDF1);
      instance.SetNewArray(&newArray_PDF1);
      instance.SetDelete(&delete_PDF1);
      instance.SetDeleteArray(&deleteArray_PDF1);
      instance.SetDestructor(&destruct_PDF1);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::PDF1*)
   {
      return GenerateInitInstanceLocal((::PDF1*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::PDF1*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_PDF2(void *p = 0);
   static void *newArray_PDF2(Long_t size, void *p);
   static void delete_PDF2(void *p);
   static void deleteArray_PDF2(void *p);
   static void destruct_PDF2(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::PDF2*)
   {
      ::PDF2 *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::PDF2 >(0);
      static ::ROOT::TGenericClassInfo 
         instance("PDF2", ::PDF2::Class_Version(), "PDF2.h", 14,
                  typeid(::PDF2), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::PDF2::Dictionary, isa_proxy, 4,
                  sizeof(::PDF2) );
      instance.SetNew(&new_PDF2);
      instance.SetNewArray(&newArray_PDF2);
      instance.SetDelete(&delete_PDF2);
      instance.SetDeleteArray(&deleteArray_PDF2);
      instance.SetDestructor(&destruct_PDF2);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::PDF2*)
   {
      return GenerateInitInstanceLocal((::PDF2*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::PDF2*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_PDF2DB(void *p = 0);
   static void *newArray_PDF2DB(Long_t size, void *p);
   static void delete_PDF2DB(void *p);
   static void deleteArray_PDF2DB(void *p);
   static void destruct_PDF2DB(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::PDF2DB*)
   {
      ::PDF2DB *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::PDF2DB >(0);
      static ::ROOT::TGenericClassInfo 
         instance("PDF2DB", ::PDF2DB::Class_Version(), "PDF2DB.h", 14,
                  typeid(::PDF2DB), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::PDF2DB::Dictionary, isa_proxy, 4,
                  sizeof(::PDF2DB) );
      instance.SetNew(&new_PDF2DB);
      instance.SetNewArray(&newArray_PDF2DB);
      instance.SetDelete(&delete_PDF2DB);
      instance.SetDeleteArray(&deleteArray_PDF2DB);
      instance.SetDestructor(&destruct_PDF2DB);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::PDF2DB*)
   {
      return GenerateInitInstanceLocal((::PDF2DB*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::PDF2DB*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_LikelihoodVar(void *p = 0);
   static void *newArray_LikelihoodVar(Long_t size, void *p);
   static void delete_LikelihoodVar(void *p);
   static void deleteArray_LikelihoodVar(void *p);
   static void destruct_LikelihoodVar(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::LikelihoodVar*)
   {
      ::LikelihoodVar *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::LikelihoodVar >(0);
      static ::ROOT::TGenericClassInfo 
         instance("LikelihoodVar", ::LikelihoodVar::Class_Version(), "LikelihoodVar.h", 14,
                  typeid(::LikelihoodVar), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::LikelihoodVar::Dictionary, isa_proxy, 4,
                  sizeof(::LikelihoodVar) );
      instance.SetNew(&new_LikelihoodVar);
      instance.SetNewArray(&newArray_LikelihoodVar);
      instance.SetDelete(&delete_LikelihoodVar);
      instance.SetDeleteArray(&deleteArray_LikelihoodVar);
      instance.SetDestructor(&destruct_LikelihoodVar);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::LikelihoodVar*)
   {
      return GenerateInitInstanceLocal((::LikelihoodVar*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::LikelihoodVar*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void delete_Likelihood(void *p);
   static void deleteArray_Likelihood(void *p);
   static void destruct_Likelihood(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::Likelihood*)
   {
      ::Likelihood *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::Likelihood >(0);
      static ::ROOT::TGenericClassInfo 
         instance("Likelihood", ::Likelihood::Class_Version(), "Likelihood.h", 16,
                  typeid(::Likelihood), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::Likelihood::Dictionary, isa_proxy, 4,
                  sizeof(::Likelihood) );
      instance.SetDelete(&delete_Likelihood);
      instance.SetDeleteArray(&deleteArray_Likelihood);
      instance.SetDestructor(&destruct_Likelihood);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::Likelihood*)
   {
      return GenerateInitInstanceLocal((::Likelihood*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::Likelihood*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static TClass *ClassifierData_Dictionary();
   static void ClassifierData_TClassManip(TClass*);
   static void *new_ClassifierData(void *p = 0);
   static void *newArray_ClassifierData(Long_t size, void *p);
   static void delete_ClassifierData(void *p);
   static void deleteArray_ClassifierData(void *p);
   static void destruct_ClassifierData(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::ClassifierData*)
   {
      ::ClassifierData *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::ClassifierData));
      static ::ROOT::TGenericClassInfo 
         instance("ClassifierData", "Classifier.h", 96,
                  typeid(::ClassifierData), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &ClassifierData_Dictionary, isa_proxy, 4,
                  sizeof(::ClassifierData) );
      instance.SetNew(&new_ClassifierData);
      instance.SetNewArray(&newArray_ClassifierData);
      instance.SetDelete(&delete_ClassifierData);
      instance.SetDeleteArray(&deleteArray_ClassifierData);
      instance.SetDestructor(&destruct_ClassifierData);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::ClassifierData*)
   {
      return GenerateInitInstanceLocal((::ClassifierData*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::ClassifierData*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *ClassifierData_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::ClassifierData*)0x0)->GetClass();
      ClassifierData_TClassManip(theClass);
   return theClass;
   }

   static void ClassifierData_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *ClassifierManager_Dictionary();
   static void ClassifierManager_TClassManip(TClass*);
   static void *new_ClassifierManager(void *p = 0);
   static void *newArray_ClassifierManager(Long_t size, void *p);
   static void delete_ClassifierManager(void *p);
   static void deleteArray_ClassifierManager(void *p);
   static void destruct_ClassifierManager(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::ClassifierManager*)
   {
      ::ClassifierManager *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::ClassifierManager));
      static ::ROOT::TGenericClassInfo 
         instance("ClassifierManager", "Classifier.h", 180,
                  typeid(::ClassifierManager), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &ClassifierManager_Dictionary, isa_proxy, 4,
                  sizeof(::ClassifierManager) );
      instance.SetNew(&new_ClassifierManager);
      instance.SetNewArray(&newArray_ClassifierManager);
      instance.SetDelete(&delete_ClassifierManager);
      instance.SetDeleteArray(&deleteArray_ClassifierManager);
      instance.SetDestructor(&destruct_ClassifierManager);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::ClassifierManager*)
   {
      return GenerateInitInstanceLocal((::ClassifierManager*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::ClassifierManager*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *ClassifierManager_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::ClassifierManager*)0x0)->GetClass();
      ClassifierManager_TClassManip(theClass);
   return theClass;
   }

   static void ClassifierManager_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *RichBDTData_Dictionary();
   static void RichBDTData_TClassManip(TClass*);
   static void *new_RichBDTData(void *p = 0);
   static void *newArray_RichBDTData(Long_t size, void *p);
   static void delete_RichBDTData(void *p);
   static void deleteArray_RichBDTData(void *p);
   static void destruct_RichBDTData(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::RichBDTData*)
   {
      ::RichBDTData *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::RichBDTData));
      static ::ROOT::TGenericClassInfo 
         instance("RichBDTData", "RichBDT.h", 15,
                  typeid(::RichBDTData), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &RichBDTData_Dictionary, isa_proxy, 4,
                  sizeof(::RichBDTData) );
      instance.SetNew(&new_RichBDTData);
      instance.SetNewArray(&newArray_RichBDTData);
      instance.SetDelete(&delete_RichBDTData);
      instance.SetDeleteArray(&deleteArray_RichBDTData);
      instance.SetDestructor(&destruct_RichBDTData);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::RichBDTData*)
   {
      return GenerateInitInstanceLocal((::RichBDTData*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::RichBDTData*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *RichBDTData_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::RichBDTData*)0x0)->GetClass();
      RichBDTData_TClassManip(theClass);
   return theClass;
   }

   static void RichBDTData_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *RichBDTMgr_Dictionary();
   static void RichBDTMgr_TClassManip(TClass*);
   static void *new_RichBDTMgr(void *p = 0);
   static void *newArray_RichBDTMgr(Long_t size, void *p);
   static void delete_RichBDTMgr(void *p);
   static void deleteArray_RichBDTMgr(void *p);
   static void destruct_RichBDTMgr(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::RichBDTMgr*)
   {
      ::RichBDTMgr *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::RichBDTMgr));
      static ::ROOT::TGenericClassInfo 
         instance("RichBDTMgr", "RichBDT.h", 52,
                  typeid(::RichBDTMgr), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &RichBDTMgr_Dictionary, isa_proxy, 4,
                  sizeof(::RichBDTMgr) );
      instance.SetNew(&new_RichBDTMgr);
      instance.SetNewArray(&newArray_RichBDTMgr);
      instance.SetDelete(&delete_RichBDTMgr);
      instance.SetDeleteArray(&deleteArray_RichBDTMgr);
      instance.SetDestructor(&destruct_RichBDTMgr);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::RichBDTMgr*)
   {
      return GenerateInitInstanceLocal((::RichBDTMgr*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::RichBDTMgr*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *RichBDTMgr_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::RichBDTMgr*)0x0)->GetClass();
      RichBDTMgr_TClassManip(theClass);
   return theClass;
   }

   static void RichBDTMgr_TClassManip(TClass* ){
   }

} // end of namespace ROOT

//______________________________________________________________________________
atomic_TClass_ptr RTIInfo::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *RTIInfo::Class_Name()
{
   return "RTIInfo";
}

//______________________________________________________________________________
const char *RTIInfo::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::RTIInfo*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int RTIInfo::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::RTIInfo*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *RTIInfo::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::RTIInfo*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *RTIInfo::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::RTIInfo*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr FileInfo::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *FileInfo::Class_Name()
{
   return "FileInfo";
}

//______________________________________________________________________________
const char *FileInfo::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::FileInfo*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int FileInfo::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::FileInfo*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *FileInfo::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::FileInfo*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *FileInfo::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::FileInfo*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr FileMCInfo::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *FileMCInfo::Class_Name()
{
   return "FileMCInfo";
}

//______________________________________________________________________________
const char *FileMCInfo::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::FileMCInfo*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int FileMCInfo::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::FileMCInfo*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *FileMCInfo::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::FileMCInfo*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *FileMCInfo::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::FileMCInfo*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr ProcInfo::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *ProcInfo::Class_Name()
{
   return "ProcInfo";
}

//______________________________________________________________________________
const char *ProcInfo::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::ProcInfo*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int ProcInfo::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::ProcInfo*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *ProcInfo::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::ProcInfo*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *ProcInfo::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::ProcInfo*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr NtpSHeader::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *NtpSHeader::Class_Name()
{
   return "NtpSHeader";
}

//______________________________________________________________________________
const char *NtpSHeader::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::NtpSHeader*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int NtpSHeader::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::NtpSHeader*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *NtpSHeader::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::NtpSHeader*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *NtpSHeader::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::NtpSHeader*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr NtpHeader::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *NtpHeader::Class_Name()
{
   return "NtpHeader";
}

//______________________________________________________________________________
const char *NtpHeader::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::NtpHeader*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int NtpHeader::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::NtpHeader*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *NtpHeader::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::NtpHeader*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *NtpHeader::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::NtpHeader*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr NtpMCHeader::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *NtpMCHeader::Class_Name()
{
   return "NtpMCHeader";
}

//______________________________________________________________________________
const char *NtpMCHeader::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::NtpMCHeader*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int NtpMCHeader::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::NtpMCHeader*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *NtpMCHeader::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::NtpMCHeader*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *NtpMCHeader::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::NtpMCHeader*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr NtpTrd::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *NtpTrd::Class_Name()
{
   return "NtpTrd";
}

//______________________________________________________________________________
const char *NtpTrd::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::NtpTrd*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int NtpTrd::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::NtpTrd*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *NtpTrd::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::NtpTrd*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *NtpTrd::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::NtpTrd*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr NtpTof::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *NtpTof::Class_Name()
{
   return "NtpTof";
}

//______________________________________________________________________________
const char *NtpTof::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::NtpTof*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int NtpTof::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::NtpTof*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *NtpTof::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::NtpTof*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *NtpTof::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::NtpTof*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr NtpTracker::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *NtpTracker::Class_Name()
{
   return "NtpTracker";
}

//______________________________________________________________________________
const char *NtpTracker::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::NtpTracker*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int NtpTracker::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::NtpTracker*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *NtpTracker::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::NtpTracker*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *NtpTracker::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::NtpTracker*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr NtpRich::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *NtpRich::Class_Name()
{
   return "NtpRich";
}

//______________________________________________________________________________
const char *NtpRich::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::NtpRich*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int NtpRich::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::NtpRich*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *NtpRich::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::NtpRich*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *NtpRich::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::NtpRich*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr NtpEcal::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *NtpEcal::Class_Name()
{
   return "NtpEcal";
}

//______________________________________________________________________________
const char *NtpEcal::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::NtpEcal*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int NtpEcal::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::NtpEcal*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *NtpEcal::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::NtpEcal*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *NtpEcal::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::NtpEcal*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr NtpAnti::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *NtpAnti::Class_Name()
{
   return "NtpAnti";
}

//______________________________________________________________________________
const char *NtpAnti::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::NtpAnti*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int NtpAnti::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::NtpAnti*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *NtpAnti::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::NtpAnti*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *NtpAnti::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::NtpAnti*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr NtpStandAlone::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *NtpStandAlone::Class_Name()
{
   return "NtpStandAlone";
}

//______________________________________________________________________________
const char *NtpStandAlone::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::NtpStandAlone*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int NtpStandAlone::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::NtpStandAlone*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *NtpStandAlone::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::NtpStandAlone*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *NtpStandAlone::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::NtpStandAlone*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr Event::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *Event::Class_Name()
{
   return "Event";
}

//______________________________________________________________________________
const char *Event::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::Event*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int Event::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::Event*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *Event::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::Event*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *Event::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::Event*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr NtpCompact::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *NtpCompact::Class_Name()
{
   return "NtpCompact";
}

//______________________________________________________________________________
const char *NtpCompact::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::NtpCompact*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int NtpCompact::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::NtpCompact*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *NtpCompact::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::NtpCompact*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *NtpCompact::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::NtpCompact*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr VPSCategory::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *VPSCategory::Class_Name()
{
   return "VPSCategory";
}

//______________________________________________________________________________
const char *VPSCategory::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::VPSCategory*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int VPSCategory::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::VPSCategory*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *VPSCategory::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::VPSCategory*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *VPSCategory::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::VPSCategory*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr VPSArchive::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *VPSArchive::Class_Name()
{
   return "VPSArchive";
}

//______________________________________________________________________________
const char *VPSArchive::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::VPSArchive*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int VPSArchive::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::VPSArchive*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *VPSArchive::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::VPSArchive*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *VPSArchive::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::VPSArchive*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr RichOccupancy::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *RichOccupancy::Class_Name()
{
   return "RichOccupancy";
}

//______________________________________________________________________________
const char *RichOccupancy::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::RichOccupancy*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int RichOccupancy::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::RichOccupancy*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *RichOccupancy::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::RichOccupancy*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *RichOccupancy::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::RichOccupancy*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr IsoPoly::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *IsoPoly::Class_Name()
{
   return "IsoPoly";
}

//______________________________________________________________________________
const char *IsoPoly::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::IsoPoly*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int IsoPoly::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::IsoPoly*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *IsoPoly::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::IsoPoly*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *IsoPoly::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::IsoPoly*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr MSpline::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *MSpline::Class_Name()
{
   return "MSpline";
}

//______________________________________________________________________________
const char *MSpline::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::MSpline*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int MSpline::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::MSpline*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *MSpline::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::MSpline*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *MSpline::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::MSpline*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr PDF1::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *PDF1::Class_Name()
{
   return "PDF1";
}

//______________________________________________________________________________
const char *PDF1::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::PDF1*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int PDF1::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::PDF1*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *PDF1::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::PDF1*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *PDF1::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::PDF1*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr PDF2::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *PDF2::Class_Name()
{
   return "PDF2";
}

//______________________________________________________________________________
const char *PDF2::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::PDF2*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int PDF2::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::PDF2*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *PDF2::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::PDF2*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *PDF2::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::PDF2*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr PDF2DB::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *PDF2DB::Class_Name()
{
   return "PDF2DB";
}

//______________________________________________________________________________
const char *PDF2DB::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::PDF2DB*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int PDF2DB::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::PDF2DB*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *PDF2DB::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::PDF2DB*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *PDF2DB::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::PDF2DB*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr LikelihoodVar::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *LikelihoodVar::Class_Name()
{
   return "LikelihoodVar";
}

//______________________________________________________________________________
const char *LikelihoodVar::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::LikelihoodVar*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int LikelihoodVar::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::LikelihoodVar*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *LikelihoodVar::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::LikelihoodVar*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *LikelihoodVar::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::LikelihoodVar*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr Likelihood::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *Likelihood::Class_Name()
{
   return "Likelihood";
}

//______________________________________________________________________________
const char *Likelihood::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::Likelihood*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int Likelihood::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::Likelihood*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *Likelihood::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::Likelihood*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *Likelihood::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::Likelihood*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
void RTIInfo::Streamer(TBuffer &R__b)
{
   // Stream an object of class RTIInfo.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(RTIInfo::Class(),this);
   } else {
      R__b.WriteClassBuffer(RTIInfo::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_RTIInfo(void *p) {
      return  p ? new(p) ::RTIInfo : new ::RTIInfo;
   }
   static void *newArray_RTIInfo(Long_t nElements, void *p) {
      return p ? new(p) ::RTIInfo[nElements] : new ::RTIInfo[nElements];
   }
   // Wrapper around operator delete
   static void delete_RTIInfo(void *p) {
      delete ((::RTIInfo*)p);
   }
   static void deleteArray_RTIInfo(void *p) {
      delete [] ((::RTIInfo*)p);
   }
   static void destruct_RTIInfo(void *p) {
      typedef ::RTIInfo current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::RTIInfo

//______________________________________________________________________________
void FileInfo::Streamer(TBuffer &R__b)
{
   // Stream an object of class FileInfo.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(FileInfo::Class(),this);
   } else {
      R__b.WriteClassBuffer(FileInfo::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_FileInfo(void *p) {
      return  p ? new(p) ::FileInfo : new ::FileInfo;
   }
   static void *newArray_FileInfo(Long_t nElements, void *p) {
      return p ? new(p) ::FileInfo[nElements] : new ::FileInfo[nElements];
   }
   // Wrapper around operator delete
   static void delete_FileInfo(void *p) {
      delete ((::FileInfo*)p);
   }
   static void deleteArray_FileInfo(void *p) {
      delete [] ((::FileInfo*)p);
   }
   static void destruct_FileInfo(void *p) {
      typedef ::FileInfo current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::FileInfo

//______________________________________________________________________________
void FileMCInfo::Streamer(TBuffer &R__b)
{
   // Stream an object of class FileMCInfo.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(FileMCInfo::Class(),this);
   } else {
      R__b.WriteClassBuffer(FileMCInfo::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_FileMCInfo(void *p) {
      return  p ? new(p) ::FileMCInfo : new ::FileMCInfo;
   }
   static void *newArray_FileMCInfo(Long_t nElements, void *p) {
      return p ? new(p) ::FileMCInfo[nElements] : new ::FileMCInfo[nElements];
   }
   // Wrapper around operator delete
   static void delete_FileMCInfo(void *p) {
      delete ((::FileMCInfo*)p);
   }
   static void deleteArray_FileMCInfo(void *p) {
      delete [] ((::FileMCInfo*)p);
   }
   static void destruct_FileMCInfo(void *p) {
      typedef ::FileMCInfo current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::FileMCInfo

//______________________________________________________________________________
void ProcInfo::Streamer(TBuffer &R__b)
{
   // Stream an object of class ProcInfo.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(ProcInfo::Class(),this);
   } else {
      R__b.WriteClassBuffer(ProcInfo::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_ProcInfo(void *p) {
      return  p ? new(p) ::ProcInfo : new ::ProcInfo;
   }
   static void *newArray_ProcInfo(Long_t nElements, void *p) {
      return p ? new(p) ::ProcInfo[nElements] : new ::ProcInfo[nElements];
   }
   // Wrapper around operator delete
   static void delete_ProcInfo(void *p) {
      delete ((::ProcInfo*)p);
   }
   static void deleteArray_ProcInfo(void *p) {
      delete [] ((::ProcInfo*)p);
   }
   static void destruct_ProcInfo(void *p) {
      typedef ::ProcInfo current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::ProcInfo

//______________________________________________________________________________
void NtpSHeader::Streamer(TBuffer &R__b)
{
   // Stream an object of class NtpSHeader.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(NtpSHeader::Class(),this);
   } else {
      R__b.WriteClassBuffer(NtpSHeader::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_NtpSHeader(void *p) {
      return  p ? new(p) ::NtpSHeader : new ::NtpSHeader;
   }
   static void *newArray_NtpSHeader(Long_t nElements, void *p) {
      return p ? new(p) ::NtpSHeader[nElements] : new ::NtpSHeader[nElements];
   }
   // Wrapper around operator delete
   static void delete_NtpSHeader(void *p) {
      delete ((::NtpSHeader*)p);
   }
   static void deleteArray_NtpSHeader(void *p) {
      delete [] ((::NtpSHeader*)p);
   }
   static void destruct_NtpSHeader(void *p) {
      typedef ::NtpSHeader current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::NtpSHeader

//______________________________________________________________________________
void NtpHeader::Streamer(TBuffer &R__b)
{
   // Stream an object of class NtpHeader.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(NtpHeader::Class(),this);
   } else {
      R__b.WriteClassBuffer(NtpHeader::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_NtpHeader(void *p) {
      return  p ? new(p) ::NtpHeader : new ::NtpHeader;
   }
   static void *newArray_NtpHeader(Long_t nElements, void *p) {
      return p ? new(p) ::NtpHeader[nElements] : new ::NtpHeader[nElements];
   }
   // Wrapper around operator delete
   static void delete_NtpHeader(void *p) {
      delete ((::NtpHeader*)p);
   }
   static void deleteArray_NtpHeader(void *p) {
      delete [] ((::NtpHeader*)p);
   }
   static void destruct_NtpHeader(void *p) {
      typedef ::NtpHeader current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::NtpHeader

//______________________________________________________________________________
void NtpMCHeader::Streamer(TBuffer &R__b)
{
   // Stream an object of class NtpMCHeader.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(NtpMCHeader::Class(),this);
   } else {
      R__b.WriteClassBuffer(NtpMCHeader::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_NtpMCHeader(void *p) {
      return  p ? new(p) ::NtpMCHeader : new ::NtpMCHeader;
   }
   static void *newArray_NtpMCHeader(Long_t nElements, void *p) {
      return p ? new(p) ::NtpMCHeader[nElements] : new ::NtpMCHeader[nElements];
   }
   // Wrapper around operator delete
   static void delete_NtpMCHeader(void *p) {
      delete ((::NtpMCHeader*)p);
   }
   static void deleteArray_NtpMCHeader(void *p) {
      delete [] ((::NtpMCHeader*)p);
   }
   static void destruct_NtpMCHeader(void *p) {
      typedef ::NtpMCHeader current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::NtpMCHeader

//______________________________________________________________________________
void NtpTrd::Streamer(TBuffer &R__b)
{
   // Stream an object of class NtpTrd.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(NtpTrd::Class(),this);
   } else {
      R__b.WriteClassBuffer(NtpTrd::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_NtpTrd(void *p) {
      return  p ? new(p) ::NtpTrd : new ::NtpTrd;
   }
   static void *newArray_NtpTrd(Long_t nElements, void *p) {
      return p ? new(p) ::NtpTrd[nElements] : new ::NtpTrd[nElements];
   }
   // Wrapper around operator delete
   static void delete_NtpTrd(void *p) {
      delete ((::NtpTrd*)p);
   }
   static void deleteArray_NtpTrd(void *p) {
      delete [] ((::NtpTrd*)p);
   }
   static void destruct_NtpTrd(void *p) {
      typedef ::NtpTrd current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::NtpTrd

//______________________________________________________________________________
void NtpTof::Streamer(TBuffer &R__b)
{
   // Stream an object of class NtpTof.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(NtpTof::Class(),this);
   } else {
      R__b.WriteClassBuffer(NtpTof::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_NtpTof(void *p) {
      return  p ? new(p) ::NtpTof : new ::NtpTof;
   }
   static void *newArray_NtpTof(Long_t nElements, void *p) {
      return p ? new(p) ::NtpTof[nElements] : new ::NtpTof[nElements];
   }
   // Wrapper around operator delete
   static void delete_NtpTof(void *p) {
      delete ((::NtpTof*)p);
   }
   static void deleteArray_NtpTof(void *p) {
      delete [] ((::NtpTof*)p);
   }
   static void destruct_NtpTof(void *p) {
      typedef ::NtpTof current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::NtpTof

//______________________________________________________________________________
void NtpTracker::Streamer(TBuffer &R__b)
{
   // Stream an object of class NtpTracker.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(NtpTracker::Class(),this);
   } else {
      R__b.WriteClassBuffer(NtpTracker::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_NtpTracker(void *p) {
      return  p ? new(p) ::NtpTracker : new ::NtpTracker;
   }
   static void *newArray_NtpTracker(Long_t nElements, void *p) {
      return p ? new(p) ::NtpTracker[nElements] : new ::NtpTracker[nElements];
   }
   // Wrapper around operator delete
   static void delete_NtpTracker(void *p) {
      delete ((::NtpTracker*)p);
   }
   static void deleteArray_NtpTracker(void *p) {
      delete [] ((::NtpTracker*)p);
   }
   static void destruct_NtpTracker(void *p) {
      typedef ::NtpTracker current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::NtpTracker

//______________________________________________________________________________
void NtpRich::Streamer(TBuffer &R__b)
{
   // Stream an object of class NtpRich.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(NtpRich::Class(),this);
   } else {
      R__b.WriteClassBuffer(NtpRich::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_NtpRich(void *p) {
      return  p ? new(p) ::NtpRich : new ::NtpRich;
   }
   static void *newArray_NtpRich(Long_t nElements, void *p) {
      return p ? new(p) ::NtpRich[nElements] : new ::NtpRich[nElements];
   }
   // Wrapper around operator delete
   static void delete_NtpRich(void *p) {
      delete ((::NtpRich*)p);
   }
   static void deleteArray_NtpRich(void *p) {
      delete [] ((::NtpRich*)p);
   }
   static void destruct_NtpRich(void *p) {
      typedef ::NtpRich current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::NtpRich

//______________________________________________________________________________
void NtpEcal::Streamer(TBuffer &R__b)
{
   // Stream an object of class NtpEcal.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(NtpEcal::Class(),this);
   } else {
      R__b.WriteClassBuffer(NtpEcal::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_NtpEcal(void *p) {
      return  p ? new(p) ::NtpEcal : new ::NtpEcal;
   }
   static void *newArray_NtpEcal(Long_t nElements, void *p) {
      return p ? new(p) ::NtpEcal[nElements] : new ::NtpEcal[nElements];
   }
   // Wrapper around operator delete
   static void delete_NtpEcal(void *p) {
      delete ((::NtpEcal*)p);
   }
   static void deleteArray_NtpEcal(void *p) {
      delete [] ((::NtpEcal*)p);
   }
   static void destruct_NtpEcal(void *p) {
      typedef ::NtpEcal current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::NtpEcal

//______________________________________________________________________________
void NtpAnti::Streamer(TBuffer &R__b)
{
   // Stream an object of class NtpAnti.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(NtpAnti::Class(),this);
   } else {
      R__b.WriteClassBuffer(NtpAnti::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_NtpAnti(void *p) {
      return  p ? new(p) ::NtpAnti : new ::NtpAnti;
   }
   static void *newArray_NtpAnti(Long_t nElements, void *p) {
      return p ? new(p) ::NtpAnti[nElements] : new ::NtpAnti[nElements];
   }
   // Wrapper around operator delete
   static void delete_NtpAnti(void *p) {
      delete ((::NtpAnti*)p);
   }
   static void deleteArray_NtpAnti(void *p) {
      delete [] ((::NtpAnti*)p);
   }
   static void destruct_NtpAnti(void *p) {
      typedef ::NtpAnti current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::NtpAnti

//______________________________________________________________________________
void NtpStandAlone::Streamer(TBuffer &R__b)
{
   // Stream an object of class NtpStandAlone.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(NtpStandAlone::Class(),this);
   } else {
      R__b.WriteClassBuffer(NtpStandAlone::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_NtpStandAlone(void *p) {
      return  p ? new(p) ::NtpStandAlone : new ::NtpStandAlone;
   }
   static void *newArray_NtpStandAlone(Long_t nElements, void *p) {
      return p ? new(p) ::NtpStandAlone[nElements] : new ::NtpStandAlone[nElements];
   }
   // Wrapper around operator delete
   static void delete_NtpStandAlone(void *p) {
      delete ((::NtpStandAlone*)p);
   }
   static void deleteArray_NtpStandAlone(void *p) {
      delete [] ((::NtpStandAlone*)p);
   }
   static void destruct_NtpStandAlone(void *p) {
      typedef ::NtpStandAlone current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::NtpStandAlone

//______________________________________________________________________________
void Event::Streamer(TBuffer &R__b)
{
   // Stream an object of class Event.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(Event::Class(),this);
   } else {
      R__b.WriteClassBuffer(Event::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_Event(void *p) {
      return  p ? new(p) ::Event : new ::Event;
   }
   static void *newArray_Event(Long_t nElements, void *p) {
      return p ? new(p) ::Event[nElements] : new ::Event[nElements];
   }
   // Wrapper around operator delete
   static void delete_Event(void *p) {
      delete ((::Event*)p);
   }
   static void deleteArray_Event(void *p) {
      delete [] ((::Event*)p);
   }
   static void destruct_Event(void *p) {
      typedef ::Event current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::Event

//______________________________________________________________________________
void NtpCompact::Streamer(TBuffer &R__b)
{
   // Stream an object of class NtpCompact.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(NtpCompact::Class(),this);
   } else {
      R__b.WriteClassBuffer(NtpCompact::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_NtpCompact(void *p) {
      return  p ? new(p) ::NtpCompact : new ::NtpCompact;
   }
   static void *newArray_NtpCompact(Long_t nElements, void *p) {
      return p ? new(p) ::NtpCompact[nElements] : new ::NtpCompact[nElements];
   }
   // Wrapper around operator delete
   static void delete_NtpCompact(void *p) {
      delete ((::NtpCompact*)p);
   }
   static void deleteArray_NtpCompact(void *p) {
      delete [] ((::NtpCompact*)p);
   }
   static void destruct_NtpCompact(void *p) {
      typedef ::NtpCompact current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::NtpCompact

//______________________________________________________________________________
void VPSCategory::Streamer(TBuffer &R__b)
{
   // Stream an object of class VPSCategory.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(VPSCategory::Class(),this);
   } else {
      R__b.WriteClassBuffer(VPSCategory::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_VPSCategory(void *p) {
      return  p ? new(p) ::VPSCategory : new ::VPSCategory;
   }
   static void *newArray_VPSCategory(Long_t nElements, void *p) {
      return p ? new(p) ::VPSCategory[nElements] : new ::VPSCategory[nElements];
   }
   // Wrapper around operator delete
   static void delete_VPSCategory(void *p) {
      delete ((::VPSCategory*)p);
   }
   static void deleteArray_VPSCategory(void *p) {
      delete [] ((::VPSCategory*)p);
   }
   static void destruct_VPSCategory(void *p) {
      typedef ::VPSCategory current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::VPSCategory

//______________________________________________________________________________
void VPSArchive::Streamer(TBuffer &R__b)
{
   // Stream an object of class VPSArchive.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(VPSArchive::Class(),this);
   } else {
      R__b.WriteClassBuffer(VPSArchive::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_VPSArchive(void *p) {
      return  p ? new(p) ::VPSArchive : new ::VPSArchive;
   }
   static void *newArray_VPSArchive(Long_t nElements, void *p) {
      return p ? new(p) ::VPSArchive[nElements] : new ::VPSArchive[nElements];
   }
   // Wrapper around operator delete
   static void delete_VPSArchive(void *p) {
      delete ((::VPSArchive*)p);
   }
   static void deleteArray_VPSArchive(void *p) {
      delete [] ((::VPSArchive*)p);
   }
   static void destruct_VPSArchive(void *p) {
      typedef ::VPSArchive current_t;
      ((current_t*)p)->~current_t();
   }
   // Wrapper around the merge function.
   static Long64_t  merge_VPSArchive(void *obj,TCollection *coll,TFileMergeInfo *) {
      return ((::VPSArchive*)obj)->Merge(coll);
   }
} // end of namespace ROOT for class ::VPSArchive

//______________________________________________________________________________
void RichOccupancy::Streamer(TBuffer &R__b)
{
   // Stream an object of class RichOccupancy.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(RichOccupancy::Class(),this);
   } else {
      R__b.WriteClassBuffer(RichOccupancy::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_RichOccupancy(void *p) {
      return  p ? new(p) ::RichOccupancy : new ::RichOccupancy;
   }
   static void *newArray_RichOccupancy(Long_t nElements, void *p) {
      return p ? new(p) ::RichOccupancy[nElements] : new ::RichOccupancy[nElements];
   }
   // Wrapper around operator delete
   static void delete_RichOccupancy(void *p) {
      delete ((::RichOccupancy*)p);
   }
   static void deleteArray_RichOccupancy(void *p) {
      delete [] ((::RichOccupancy*)p);
   }
   static void destruct_RichOccupancy(void *p) {
      typedef ::RichOccupancy current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::RichOccupancy

//______________________________________________________________________________
void IsoPoly::Streamer(TBuffer &R__b)
{
   // Stream an object of class IsoPoly.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(IsoPoly::Class(),this);
   } else {
      R__b.WriteClassBuffer(IsoPoly::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_IsoPoly(void *p) {
      return  p ? new(p) ::IsoPoly : new ::IsoPoly;
   }
   static void *newArray_IsoPoly(Long_t nElements, void *p) {
      return p ? new(p) ::IsoPoly[nElements] : new ::IsoPoly[nElements];
   }
   // Wrapper around operator delete
   static void delete_IsoPoly(void *p) {
      delete ((::IsoPoly*)p);
   }
   static void deleteArray_IsoPoly(void *p) {
      delete [] ((::IsoPoly*)p);
   }
   static void destruct_IsoPoly(void *p) {
      typedef ::IsoPoly current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::IsoPoly

//______________________________________________________________________________
void MSpline::Streamer(TBuffer &R__b)
{
   // Stream an object of class MSpline.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(MSpline::Class(),this);
   } else {
      R__b.WriteClassBuffer(MSpline::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_MSpline(void *p) {
      return  p ? new(p) ::MSpline : new ::MSpline;
   }
   static void *newArray_MSpline(Long_t nElements, void *p) {
      return p ? new(p) ::MSpline[nElements] : new ::MSpline[nElements];
   }
   // Wrapper around operator delete
   static void delete_MSpline(void *p) {
      delete ((::MSpline*)p);
   }
   static void deleteArray_MSpline(void *p) {
      delete [] ((::MSpline*)p);
   }
   static void destruct_MSpline(void *p) {
      typedef ::MSpline current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::MSpline

//______________________________________________________________________________
void PDF1::Streamer(TBuffer &R__b)
{
   // Stream an object of class PDF1.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(PDF1::Class(),this);
   } else {
      R__b.WriteClassBuffer(PDF1::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_PDF1(void *p) {
      return  p ? new(p) ::PDF1 : new ::PDF1;
   }
   static void *newArray_PDF1(Long_t nElements, void *p) {
      return p ? new(p) ::PDF1[nElements] : new ::PDF1[nElements];
   }
   // Wrapper around operator delete
   static void delete_PDF1(void *p) {
      delete ((::PDF1*)p);
   }
   static void deleteArray_PDF1(void *p) {
      delete [] ((::PDF1*)p);
   }
   static void destruct_PDF1(void *p) {
      typedef ::PDF1 current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::PDF1

//______________________________________________________________________________
void PDF2::Streamer(TBuffer &R__b)
{
   // Stream an object of class PDF2.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(PDF2::Class(),this);
   } else {
      R__b.WriteClassBuffer(PDF2::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_PDF2(void *p) {
      return  p ? new(p) ::PDF2 : new ::PDF2;
   }
   static void *newArray_PDF2(Long_t nElements, void *p) {
      return p ? new(p) ::PDF2[nElements] : new ::PDF2[nElements];
   }
   // Wrapper around operator delete
   static void delete_PDF2(void *p) {
      delete ((::PDF2*)p);
   }
   static void deleteArray_PDF2(void *p) {
      delete [] ((::PDF2*)p);
   }
   static void destruct_PDF2(void *p) {
      typedef ::PDF2 current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::PDF2

//______________________________________________________________________________
void PDF2DB::Streamer(TBuffer &R__b)
{
   // Stream an object of class PDF2DB.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(PDF2DB::Class(),this);
   } else {
      R__b.WriteClassBuffer(PDF2DB::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_PDF2DB(void *p) {
      return  p ? new(p) ::PDF2DB : new ::PDF2DB;
   }
   static void *newArray_PDF2DB(Long_t nElements, void *p) {
      return p ? new(p) ::PDF2DB[nElements] : new ::PDF2DB[nElements];
   }
   // Wrapper around operator delete
   static void delete_PDF2DB(void *p) {
      delete ((::PDF2DB*)p);
   }
   static void deleteArray_PDF2DB(void *p) {
      delete [] ((::PDF2DB*)p);
   }
   static void destruct_PDF2DB(void *p) {
      typedef ::PDF2DB current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::PDF2DB

//______________________________________________________________________________
void LikelihoodVar::Streamer(TBuffer &R__b)
{
   // Stream an object of class LikelihoodVar.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(LikelihoodVar::Class(),this);
   } else {
      R__b.WriteClassBuffer(LikelihoodVar::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_LikelihoodVar(void *p) {
      return  p ? new(p) ::LikelihoodVar : new ::LikelihoodVar;
   }
   static void *newArray_LikelihoodVar(Long_t nElements, void *p) {
      return p ? new(p) ::LikelihoodVar[nElements] : new ::LikelihoodVar[nElements];
   }
   // Wrapper around operator delete
   static void delete_LikelihoodVar(void *p) {
      delete ((::LikelihoodVar*)p);
   }
   static void deleteArray_LikelihoodVar(void *p) {
      delete [] ((::LikelihoodVar*)p);
   }
   static void destruct_LikelihoodVar(void *p) {
      typedef ::LikelihoodVar current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::LikelihoodVar

//______________________________________________________________________________
void Likelihood::Streamer(TBuffer &R__b)
{
   // Stream an object of class Likelihood.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(Likelihood::Class(),this);
   } else {
      R__b.WriteClassBuffer(Likelihood::Class(),this);
   }
}

namespace ROOT {
   // Wrapper around operator delete
   static void delete_Likelihood(void *p) {
      delete ((::Likelihood*)p);
   }
   static void deleteArray_Likelihood(void *p) {
      delete [] ((::Likelihood*)p);
   }
   static void destruct_Likelihood(void *p) {
      typedef ::Likelihood current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::Likelihood

namespace ROOT {
   // Wrappers around operator new
   static void *new_ClassifierData(void *p) {
      return  p ? new(p) ::ClassifierData : new ::ClassifierData;
   }
   static void *newArray_ClassifierData(Long_t nElements, void *p) {
      return p ? new(p) ::ClassifierData[nElements] : new ::ClassifierData[nElements];
   }
   // Wrapper around operator delete
   static void delete_ClassifierData(void *p) {
      delete ((::ClassifierData*)p);
   }
   static void deleteArray_ClassifierData(void *p) {
      delete [] ((::ClassifierData*)p);
   }
   static void destruct_ClassifierData(void *p) {
      typedef ::ClassifierData current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::ClassifierData

namespace ROOT {
   // Wrappers around operator new
   static void *new_ClassifierManager(void *p) {
      return  p ? new(p) ::ClassifierManager : new ::ClassifierManager;
   }
   static void *newArray_ClassifierManager(Long_t nElements, void *p) {
      return p ? new(p) ::ClassifierManager[nElements] : new ::ClassifierManager[nElements];
   }
   // Wrapper around operator delete
   static void delete_ClassifierManager(void *p) {
      delete ((::ClassifierManager*)p);
   }
   static void deleteArray_ClassifierManager(void *p) {
      delete [] ((::ClassifierManager*)p);
   }
   static void destruct_ClassifierManager(void *p) {
      typedef ::ClassifierManager current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::ClassifierManager

namespace ROOT {
   // Wrappers around operator new
   static void *new_RichBDTData(void *p) {
      return  p ? new(p) ::RichBDTData : new ::RichBDTData;
   }
   static void *newArray_RichBDTData(Long_t nElements, void *p) {
      return p ? new(p) ::RichBDTData[nElements] : new ::RichBDTData[nElements];
   }
   // Wrapper around operator delete
   static void delete_RichBDTData(void *p) {
      delete ((::RichBDTData*)p);
   }
   static void deleteArray_RichBDTData(void *p) {
      delete [] ((::RichBDTData*)p);
   }
   static void destruct_RichBDTData(void *p) {
      typedef ::RichBDTData current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::RichBDTData

namespace ROOT {
   // Wrappers around operator new
   static void *new_RichBDTMgr(void *p) {
      return  p ? new(p) ::RichBDTMgr : new ::RichBDTMgr;
   }
   static void *newArray_RichBDTMgr(Long_t nElements, void *p) {
      return p ? new(p) ::RichBDTMgr[nElements] : new ::RichBDTMgr[nElements];
   }
   // Wrapper around operator delete
   static void delete_RichBDTMgr(void *p) {
      delete ((::RichBDTMgr*)p);
   }
   static void deleteArray_RichBDTMgr(void *p) {
      delete [] ((::RichBDTMgr*)p);
   }
   static void destruct_RichBDTMgr(void *p) {
      typedef ::RichBDTMgr current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::RichBDTMgr

namespace ROOT {
   static TClass *vectorlEintgR_Dictionary();
   static void vectorlEintgR_TClassManip(TClass*);
   static void *new_vectorlEintgR(void *p = 0);
   static void *newArray_vectorlEintgR(Long_t size, void *p);
   static void delete_vectorlEintgR(void *p);
   static void deleteArray_vectorlEintgR(void *p);
   static void destruct_vectorlEintgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<int>*)
   {
      vector<int> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<int>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<int>", -2, "vector", 214,
                  typeid(vector<int>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEintgR_Dictionary, isa_proxy, 0,
                  sizeof(vector<int>) );
      instance.SetNew(&new_vectorlEintgR);
      instance.SetNewArray(&newArray_vectorlEintgR);
      instance.SetDelete(&delete_vectorlEintgR);
      instance.SetDeleteArray(&deleteArray_vectorlEintgR);
      instance.SetDestructor(&destruct_vectorlEintgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<int> >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const vector<int>*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEintgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<int>*)0x0)->GetClass();
      vectorlEintgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEintgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEintgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<int> : new vector<int>;
   }
   static void *newArray_vectorlEintgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<int>[nElements] : new vector<int>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEintgR(void *p) {
      delete ((vector<int>*)p);
   }
   static void deleteArray_vectorlEintgR(void *p) {
      delete [] ((vector<int>*)p);
   }
   static void destruct_vectorlEintgR(void *p) {
      typedef vector<int> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<int>

namespace ROOT {
   static TClass *vectorlEVPSCategorygR_Dictionary();
   static void vectorlEVPSCategorygR_TClassManip(TClass*);
   static void *new_vectorlEVPSCategorygR(void *p = 0);
   static void *newArray_vectorlEVPSCategorygR(Long_t size, void *p);
   static void delete_vectorlEVPSCategorygR(void *p);
   static void deleteArray_vectorlEVPSCategorygR(void *p);
   static void destruct_vectorlEVPSCategorygR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<VPSCategory>*)
   {
      vector<VPSCategory> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<VPSCategory>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<VPSCategory>", -2, "vector", 214,
                  typeid(vector<VPSCategory>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEVPSCategorygR_Dictionary, isa_proxy, 0,
                  sizeof(vector<VPSCategory>) );
      instance.SetNew(&new_vectorlEVPSCategorygR);
      instance.SetNewArray(&newArray_vectorlEVPSCategorygR);
      instance.SetDelete(&delete_vectorlEVPSCategorygR);
      instance.SetDeleteArray(&deleteArray_vectorlEVPSCategorygR);
      instance.SetDestructor(&destruct_vectorlEVPSCategorygR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<VPSCategory> >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const vector<VPSCategory>*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEVPSCategorygR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<VPSCategory>*)0x0)->GetClass();
      vectorlEVPSCategorygR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEVPSCategorygR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEVPSCategorygR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<VPSCategory> : new vector<VPSCategory>;
   }
   static void *newArray_vectorlEVPSCategorygR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<VPSCategory>[nElements] : new vector<VPSCategory>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEVPSCategorygR(void *p) {
      delete ((vector<VPSCategory>*)p);
   }
   static void deleteArray_vectorlEVPSCategorygR(void *p) {
      delete [] ((vector<VPSCategory>*)p);
   }
   static void destruct_vectorlEVPSCategorygR(void *p) {
      typedef vector<VPSCategory> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<VPSCategory>

namespace ROOT {
   static TClass *vectorlELikelihoodVarmUgR_Dictionary();
   static void vectorlELikelihoodVarmUgR_TClassManip(TClass*);
   static void *new_vectorlELikelihoodVarmUgR(void *p = 0);
   static void *newArray_vectorlELikelihoodVarmUgR(Long_t size, void *p);
   static void delete_vectorlELikelihoodVarmUgR(void *p);
   static void deleteArray_vectorlELikelihoodVarmUgR(void *p);
   static void destruct_vectorlELikelihoodVarmUgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<LikelihoodVar*>*)
   {
      vector<LikelihoodVar*> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<LikelihoodVar*>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<LikelihoodVar*>", -2, "vector", 214,
                  typeid(vector<LikelihoodVar*>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlELikelihoodVarmUgR_Dictionary, isa_proxy, 0,
                  sizeof(vector<LikelihoodVar*>) );
      instance.SetNew(&new_vectorlELikelihoodVarmUgR);
      instance.SetNewArray(&newArray_vectorlELikelihoodVarmUgR);
      instance.SetDelete(&delete_vectorlELikelihoodVarmUgR);
      instance.SetDeleteArray(&deleteArray_vectorlELikelihoodVarmUgR);
      instance.SetDestructor(&destruct_vectorlELikelihoodVarmUgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<LikelihoodVar*> >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const vector<LikelihoodVar*>*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlELikelihoodVarmUgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<LikelihoodVar*>*)0x0)->GetClass();
      vectorlELikelihoodVarmUgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlELikelihoodVarmUgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlELikelihoodVarmUgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<LikelihoodVar*> : new vector<LikelihoodVar*>;
   }
   static void *newArray_vectorlELikelihoodVarmUgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<LikelihoodVar*>[nElements] : new vector<LikelihoodVar*>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlELikelihoodVarmUgR(void *p) {
      delete ((vector<LikelihoodVar*>*)p);
   }
   static void deleteArray_vectorlELikelihoodVarmUgR(void *p) {
      delete [] ((vector<LikelihoodVar*>*)p);
   }
   static void destruct_vectorlELikelihoodVarmUgR(void *p) {
      typedef vector<LikelihoodVar*> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<LikelihoodVar*>

namespace ROOT {
   static TClass *maplEstringcOPDF2mUgR_Dictionary();
   static void maplEstringcOPDF2mUgR_TClassManip(TClass*);
   static void *new_maplEstringcOPDF2mUgR(void *p = 0);
   static void *newArray_maplEstringcOPDF2mUgR(Long_t size, void *p);
   static void delete_maplEstringcOPDF2mUgR(void *p);
   static void deleteArray_maplEstringcOPDF2mUgR(void *p);
   static void destruct_maplEstringcOPDF2mUgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const map<string,PDF2*>*)
   {
      map<string,PDF2*> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(map<string,PDF2*>));
      static ::ROOT::TGenericClassInfo 
         instance("map<string,PDF2*>", -2, "map", 96,
                  typeid(map<string,PDF2*>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &maplEstringcOPDF2mUgR_Dictionary, isa_proxy, 0,
                  sizeof(map<string,PDF2*>) );
      instance.SetNew(&new_maplEstringcOPDF2mUgR);
      instance.SetNewArray(&newArray_maplEstringcOPDF2mUgR);
      instance.SetDelete(&delete_maplEstringcOPDF2mUgR);
      instance.SetDeleteArray(&deleteArray_maplEstringcOPDF2mUgR);
      instance.SetDestructor(&destruct_maplEstringcOPDF2mUgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::MapInsert< map<string,PDF2*> >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const map<string,PDF2*>*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *maplEstringcOPDF2mUgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const map<string,PDF2*>*)0x0)->GetClass();
      maplEstringcOPDF2mUgR_TClassManip(theClass);
   return theClass;
   }

   static void maplEstringcOPDF2mUgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_maplEstringcOPDF2mUgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) map<string,PDF2*> : new map<string,PDF2*>;
   }
   static void *newArray_maplEstringcOPDF2mUgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) map<string,PDF2*>[nElements] : new map<string,PDF2*>[nElements];
   }
   // Wrapper around operator delete
   static void delete_maplEstringcOPDF2mUgR(void *p) {
      delete ((map<string,PDF2*>*)p);
   }
   static void deleteArray_maplEstringcOPDF2mUgR(void *p) {
      delete [] ((map<string,PDF2*>*)p);
   }
   static void destruct_maplEstringcOPDF2mUgR(void *p) {
      typedef map<string,PDF2*> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class map<string,PDF2*>

namespace ROOT {
   static TClass *maplEintcOstringgR_Dictionary();
   static void maplEintcOstringgR_TClassManip(TClass*);
   static void *new_maplEintcOstringgR(void *p = 0);
   static void *newArray_maplEintcOstringgR(Long_t size, void *p);
   static void delete_maplEintcOstringgR(void *p);
   static void deleteArray_maplEintcOstringgR(void *p);
   static void destruct_maplEintcOstringgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const map<int,string>*)
   {
      map<int,string> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(map<int,string>));
      static ::ROOT::TGenericClassInfo 
         instance("map<int,string>", -2, "map", 96,
                  typeid(map<int,string>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &maplEintcOstringgR_Dictionary, isa_proxy, 0,
                  sizeof(map<int,string>) );
      instance.SetNew(&new_maplEintcOstringgR);
      instance.SetNewArray(&newArray_maplEintcOstringgR);
      instance.SetDelete(&delete_maplEintcOstringgR);
      instance.SetDeleteArray(&deleteArray_maplEintcOstringgR);
      instance.SetDestructor(&destruct_maplEintcOstringgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::MapInsert< map<int,string> >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const map<int,string>*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *maplEintcOstringgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const map<int,string>*)0x0)->GetClass();
      maplEintcOstringgR_TClassManip(theClass);
   return theClass;
   }

   static void maplEintcOstringgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_maplEintcOstringgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) map<int,string> : new map<int,string>;
   }
   static void *newArray_maplEintcOstringgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) map<int,string>[nElements] : new map<int,string>[nElements];
   }
   // Wrapper around operator delete
   static void delete_maplEintcOstringgR(void *p) {
      delete ((map<int,string>*)p);
   }
   static void deleteArray_maplEintcOstringgR(void *p) {
      delete [] ((map<int,string>*)p);
   }
   static void destruct_maplEintcOstringgR(void *p) {
      typedef map<int,string> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class map<int,string>

namespace ROOT {
   static TClass *maplEintcOTMVAcLcLReadermUgR_Dictionary();
   static void maplEintcOTMVAcLcLReadermUgR_TClassManip(TClass*);
   static void *new_maplEintcOTMVAcLcLReadermUgR(void *p = 0);
   static void *newArray_maplEintcOTMVAcLcLReadermUgR(Long_t size, void *p);
   static void delete_maplEintcOTMVAcLcLReadermUgR(void *p);
   static void deleteArray_maplEintcOTMVAcLcLReadermUgR(void *p);
   static void destruct_maplEintcOTMVAcLcLReadermUgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const map<int,TMVA::Reader*>*)
   {
      map<int,TMVA::Reader*> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(map<int,TMVA::Reader*>));
      static ::ROOT::TGenericClassInfo 
         instance("map<int,TMVA::Reader*>", -2, "map", 96,
                  typeid(map<int,TMVA::Reader*>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &maplEintcOTMVAcLcLReadermUgR_Dictionary, isa_proxy, 0,
                  sizeof(map<int,TMVA::Reader*>) );
      instance.SetNew(&new_maplEintcOTMVAcLcLReadermUgR);
      instance.SetNewArray(&newArray_maplEintcOTMVAcLcLReadermUgR);
      instance.SetDelete(&delete_maplEintcOTMVAcLcLReadermUgR);
      instance.SetDeleteArray(&deleteArray_maplEintcOTMVAcLcLReadermUgR);
      instance.SetDestructor(&destruct_maplEintcOTMVAcLcLReadermUgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::MapInsert< map<int,TMVA::Reader*> >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const map<int,TMVA::Reader*>*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *maplEintcOTMVAcLcLReadermUgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const map<int,TMVA::Reader*>*)0x0)->GetClass();
      maplEintcOTMVAcLcLReadermUgR_TClassManip(theClass);
   return theClass;
   }

   static void maplEintcOTMVAcLcLReadermUgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_maplEintcOTMVAcLcLReadermUgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) map<int,TMVA::Reader*> : new map<int,TMVA::Reader*>;
   }
   static void *newArray_maplEintcOTMVAcLcLReadermUgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) map<int,TMVA::Reader*>[nElements] : new map<int,TMVA::Reader*>[nElements];
   }
   // Wrapper around operator delete
   static void delete_maplEintcOTMVAcLcLReadermUgR(void *p) {
      delete ((map<int,TMVA::Reader*>*)p);
   }
   static void deleteArray_maplEintcOTMVAcLcLReadermUgR(void *p) {
      delete [] ((map<int,TMVA::Reader*>*)p);
   }
   static void destruct_maplEintcOTMVAcLcLReadermUgR(void *p) {
      typedef map<int,TMVA::Reader*> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class map<int,TMVA::Reader*>

namespace ROOT {
   static TClass *bitsetlE12gR_Dictionary();
   static void bitsetlE12gR_TClassManip(TClass*);
   static void *new_bitsetlE12gR(void *p = 0);
   static void *newArray_bitsetlE12gR(Long_t size, void *p);
   static void delete_bitsetlE12gR(void *p);
   static void deleteArray_bitsetlE12gR(void *p);
   static void destruct_bitsetlE12gR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const bitset<12>*)
   {
      bitset<12> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(bitset<12>));
      static ::ROOT::TGenericClassInfo 
         instance("bitset<12>", 2, "bitset", 748,
                  typeid(bitset<12>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &bitsetlE12gR_Dictionary, isa_proxy, 0,
                  sizeof(bitset<12>) );
      instance.SetNew(&new_bitsetlE12gR);
      instance.SetNewArray(&newArray_bitsetlE12gR);
      instance.SetDelete(&delete_bitsetlE12gR);
      instance.SetDeleteArray(&deleteArray_bitsetlE12gR);
      instance.SetDestructor(&destruct_bitsetlE12gR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback<Internal::TStdBitsetHelper< bitset<12> > >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const bitset<12>*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *bitsetlE12gR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const bitset<12>*)0x0)->GetClass();
      bitsetlE12gR_TClassManip(theClass);
   return theClass;
   }

   static void bitsetlE12gR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_bitsetlE12gR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) bitset<12> : new bitset<12>;
   }
   static void *newArray_bitsetlE12gR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) bitset<12>[nElements] : new bitset<12>[nElements];
   }
   // Wrapper around operator delete
   static void delete_bitsetlE12gR(void *p) {
      delete ((bitset<12>*)p);
   }
   static void deleteArray_bitsetlE12gR(void *p) {
      delete [] ((bitset<12>*)p);
   }
   static void destruct_bitsetlE12gR(void *p) {
      typedef bitset<12> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class bitset<12>

namespace {
  void TriggerDictionaryInitialization_libntp_Impl() {
    static const char* headers[] = {
"/afs/cern.ch/user/f/fdimicco/Work/Deutons/data_reduction/include/Ntp.h",
"/afs/cern.ch/user/f/fdimicco/Work/Deutons/data_reduction/include/VPSArchive.h",
"/afs/cern.ch/user/f/fdimicco/Work/Deutons/data_reduction/include/RichOccupancy.h",
"/afs/cern.ch/user/f/fdimicco/Work/Deutons/data_reduction/include/Classifier.h",
"/afs/cern.ch/user/f/fdimicco/Work/Deutons/data_reduction/include/RichBDT.h",
"/afs/cern.ch/user/f/fdimicco/Work/Deutons/data_reduction/include/MSpline.h",
"/afs/cern.ch/user/f/fdimicco/Work/Deutons/data_reduction/include/IsoPoly.h",
"/afs/cern.ch/user/f/fdimicco/Work/Deutons/data_reduction/include/PDF1.h",
"/afs/cern.ch/user/f/fdimicco/Work/Deutons/data_reduction/include/PDF2.h",
"/afs/cern.ch/user/f/fdimicco/Work/Deutons/data_reduction/include/PDF2DB.h",
"/afs/cern.ch/user/f/fdimicco/Work/Deutons/data_reduction/include/LikelihoodVar.h",
"/afs/cern.ch/user/f/fdimicco/Work/Deutons/data_reduction/include/Likelihood.h",
0
    };
    static const char* includePaths[] = {
"/cvmfs/sft.cern.ch/lcg/app/releases/ROOT/6.14.04/x86_64-centos7-gcc48-opt/root/include",
"/afs/cern.ch/user/f/fdimicco/Work/Deutons/data_reduction/include",
"/data1/home/pzuccon/work/dbar/AMS/include",
"/cvmfs/sft.cern.ch/lcg/app/releases/ROOT/6.14.04/x86_64-centos7-gcc48-opt/root/include",
"/afs/cern.ch/work/f/fdimicco/private/Deutons/data_reduction/build/",
0
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "libntp dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_Autoloading_Map;
class __attribute__((annotate("$clingAutoload$/afs/cern.ch/user/f/fdimicco/Work/Deutons/data_reduction/include/Ntp.h")))  RTIInfo;
class __attribute__((annotate("$clingAutoload$/afs/cern.ch/user/f/fdimicco/Work/Deutons/data_reduction/include/Ntp.h")))  FileInfo;
class __attribute__((annotate("$clingAutoload$/afs/cern.ch/user/f/fdimicco/Work/Deutons/data_reduction/include/Ntp.h")))  FileMCInfo;
class __attribute__((annotate("$clingAutoload$/afs/cern.ch/user/f/fdimicco/Work/Deutons/data_reduction/include/Ntp.h")))  ProcInfo;
class __attribute__((annotate("$clingAutoload$/afs/cern.ch/user/f/fdimicco/Work/Deutons/data_reduction/include/Ntp.h")))  NtpSHeader;
class __attribute__((annotate("$clingAutoload$/afs/cern.ch/user/f/fdimicco/Work/Deutons/data_reduction/include/Ntp.h")))  NtpHeader;
class __attribute__((annotate("$clingAutoload$/afs/cern.ch/user/f/fdimicco/Work/Deutons/data_reduction/include/Ntp.h")))  NtpMCHeader;
class __attribute__((annotate("$clingAutoload$/afs/cern.ch/user/f/fdimicco/Work/Deutons/data_reduction/include/Ntp.h")))  NtpTrd;
class __attribute__((annotate("$clingAutoload$/afs/cern.ch/user/f/fdimicco/Work/Deutons/data_reduction/include/Ntp.h")))  NtpTof;
class __attribute__((annotate("$clingAutoload$/afs/cern.ch/user/f/fdimicco/Work/Deutons/data_reduction/include/Ntp.h")))  NtpTracker;
class __attribute__((annotate("$clingAutoload$/afs/cern.ch/user/f/fdimicco/Work/Deutons/data_reduction/include/Ntp.h")))  NtpRich;
class __attribute__((annotate("$clingAutoload$/afs/cern.ch/user/f/fdimicco/Work/Deutons/data_reduction/include/Ntp.h")))  NtpEcal;
class __attribute__((annotate("$clingAutoload$/afs/cern.ch/user/f/fdimicco/Work/Deutons/data_reduction/include/Ntp.h")))  NtpAnti;
class __attribute__((annotate("$clingAutoload$/afs/cern.ch/user/f/fdimicco/Work/Deutons/data_reduction/include/Ntp.h")))  NtpStandAlone;
class __attribute__((annotate("$clingAutoload$/afs/cern.ch/user/f/fdimicco/Work/Deutons/data_reduction/include/Ntp.h")))  Event;
class __attribute__((annotate("$clingAutoload$/afs/cern.ch/user/f/fdimicco/Work/Deutons/data_reduction/include/Ntp.h")))  NtpCompact;
class __attribute__((annotate("$clingAutoload$/afs/cern.ch/user/f/fdimicco/Work/Deutons/data_reduction/include/VPSArchive.h")))  VPSCategory;
class __attribute__((annotate("$clingAutoload$/afs/cern.ch/user/f/fdimicco/Work/Deutons/data_reduction/include/VPSArchive.h")))  VPSArchive;
class __attribute__((annotate("$clingAutoload$/afs/cern.ch/user/f/fdimicco/Work/Deutons/data_reduction/include/RichOccupancy.h")))  RichOccupancy;
class __attribute__((annotate("$clingAutoload$IsoPoly.h")))  __attribute__((annotate("$clingAutoload$/afs/cern.ch/user/f/fdimicco/Work/Deutons/data_reduction/include/Classifier.h")))  IsoPoly;
class __attribute__((annotate("$clingAutoload$MSpline.h")))  __attribute__((annotate("$clingAutoload$/afs/cern.ch/user/f/fdimicco/Work/Deutons/data_reduction/include/Classifier.h")))  MSpline;
class __attribute__((annotate("$clingAutoload$PDF1.h")))  __attribute__((annotate("$clingAutoload$/afs/cern.ch/user/f/fdimicco/Work/Deutons/data_reduction/include/Classifier.h")))  PDF1;
class __attribute__((annotate("$clingAutoload$PDF2.h")))  __attribute__((annotate("$clingAutoload$/afs/cern.ch/user/f/fdimicco/Work/Deutons/data_reduction/include/Classifier.h")))  PDF2;
class __attribute__((annotate("$clingAutoload$PDF2DB.h")))  __attribute__((annotate("$clingAutoload$/afs/cern.ch/user/f/fdimicco/Work/Deutons/data_reduction/include/Classifier.h")))  PDF2DB;
class __attribute__((annotate("$clingAutoload$LikelihoodVar.h")))  __attribute__((annotate("$clingAutoload$/afs/cern.ch/user/f/fdimicco/Work/Deutons/data_reduction/include/Classifier.h")))  LikelihoodVar;
class __attribute__((annotate("$clingAutoload$Likelihood.h")))  __attribute__((annotate("$clingAutoload$/afs/cern.ch/user/f/fdimicco/Work/Deutons/data_reduction/include/Classifier.h")))  Likelihood;
class __attribute__((annotate("$clingAutoload$/afs/cern.ch/user/f/fdimicco/Work/Deutons/data_reduction/include/Classifier.h")))  ClassifierData;
class __attribute__((annotate("$clingAutoload$/afs/cern.ch/user/f/fdimicco/Work/Deutons/data_reduction/include/Classifier.h")))  ClassifierManager;
class __attribute__((annotate("$clingAutoload$/afs/cern.ch/user/f/fdimicco/Work/Deutons/data_reduction/include/RichBDT.h")))  RichBDTData;
class __attribute__((annotate("$clingAutoload$/afs/cern.ch/user/f/fdimicco/Work/Deutons/data_reduction/include/RichBDT.h")))  RichBDTMgr;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "libntp dictionary payload"

#ifndef G__VECTOR_HAS_CLASS_ITERATOR
  #define G__VECTOR_HAS_CLASS_ITERATOR 1
#endif
#ifndef __ROOTSHAREDLIBRARY__
  #define __ROOTSHAREDLIBRARY__ 1
#endif
#ifndef _PGTRACK_
  #define _PGTRACK_ 1
#endif
#ifndef _NTPLIB_
  #define _NTPLIB_ 1
#endif

#define _BACKWARD_BACKWARD_WARNING_H
#include "/afs/cern.ch/user/f/fdimicco/Work/Deutons/data_reduction/include/Ntp.h"
#include "/afs/cern.ch/user/f/fdimicco/Work/Deutons/data_reduction/include/VPSArchive.h"
#include "/afs/cern.ch/user/f/fdimicco/Work/Deutons/data_reduction/include/RichOccupancy.h"
#include "/afs/cern.ch/user/f/fdimicco/Work/Deutons/data_reduction/include/Classifier.h"
#include "/afs/cern.ch/user/f/fdimicco/Work/Deutons/data_reduction/include/RichBDT.h"
#include "/afs/cern.ch/user/f/fdimicco/Work/Deutons/data_reduction/include/MSpline.h"
#include "/afs/cern.ch/user/f/fdimicco/Work/Deutons/data_reduction/include/IsoPoly.h"
#include "/afs/cern.ch/user/f/fdimicco/Work/Deutons/data_reduction/include/PDF1.h"
#include "/afs/cern.ch/user/f/fdimicco/Work/Deutons/data_reduction/include/PDF2.h"
#include "/afs/cern.ch/user/f/fdimicco/Work/Deutons/data_reduction/include/PDF2DB.h"
#include "/afs/cern.ch/user/f/fdimicco/Work/Deutons/data_reduction/include/LikelihoodVar.h"
#include "/afs/cern.ch/user/f/fdimicco/Work/Deutons/data_reduction/include/Likelihood.h"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[]={
"ClassifierData", payloadCode, "@",
"ClassifierManager", payloadCode, "@",
"Event", payloadCode, "@",
"FileInfo", payloadCode, "@",
"FileMCInfo", payloadCode, "@",
"IsoPoly", payloadCode, "@",
"Likelihood", payloadCode, "@",
"LikelihoodVar", payloadCode, "@",
"MSpline", payloadCode, "@",
"NtpAnti", payloadCode, "@",
"NtpCompact", payloadCode, "@",
"NtpEcal", payloadCode, "@",
"NtpHeader", payloadCode, "@",
"NtpMCHeader", payloadCode, "@",
"NtpRich", payloadCode, "@",
"NtpSHeader", payloadCode, "@",
"NtpStandAlone", payloadCode, "@",
"NtpTof", payloadCode, "@",
"NtpTracker", payloadCode, "@",
"NtpTrd", payloadCode, "@",
"PDF1", payloadCode, "@",
"PDF2", payloadCode, "@",
"PDF2DB", payloadCode, "@",
"ProcInfo", payloadCode, "@",
"RTIInfo", payloadCode, "@",
"RichBDTData", payloadCode, "@",
"RichBDTMgr", payloadCode, "@",
"RichOccupancy", payloadCode, "@",
"VPSArchive", payloadCode, "@",
"VPSCategory", payloadCode, "@",
nullptr};

    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("libntp",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_libntp_Impl, {}, classesHeaders, /*has no C++ module*/false);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_libntp_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_libntp() {
  TriggerDictionaryInitialization_libntp_Impl();
}
