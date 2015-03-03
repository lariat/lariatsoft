////////////////////////////////////////////////////////////////////////
// Class:       DigitReader
// Module Type: analyzer
// File:        DigitReader_module.cc
//
// Generated at Tue Feb 10 15:40:09 2015 by Will Flanagan using artmod
// from cetpkgsupport v1_08_02.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "RawData/AuxDetDigit.h"

#include <memory>
#include <iostream>
#include <vector>

class DigitReader;

class DigitReader : public art::EDAnalyzer {
public:
  explicit DigitReader(fhicl::ParameterSet const & p);
  // The destructor generated by the compiler is fine for classes
  // without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  DigitReader(DigitReader const &) = delete;
  DigitReader(DigitReader &&) = delete;
  DigitReader & operator = (DigitReader const &) = delete;
  DigitReader & operator = (DigitReader &&) = delete;

  // Required functions.
  void analyze(art::Event const & e) override;

  // Selected optional functions.
  void beginJob() override;
  void reconfigure(fhicl::ParameterSet const & p) override;

private:

  std::string fCaenV1740Board8Label;
  std::string fCaenV1751Board1Label;
  std::string fCaenV1751Board2Label;

};

//------------------------------------------------------------------------------
DigitReader::DigitReader(fhicl::ParameterSet const & p)
  : EDAnalyzer(p)
{
  this->reconfigure(p);
}

//------------------------------------------------------------------------------
void DigitReader::reconfigure(fhicl::ParameterSet const & p)
{
  fCaenV1740Board8Label = p.get< std::string >("CaenV1740Board8Label",
                                               "CaenV1740Board8");
  fCaenV1751Board1Label = p.get< std::string >("CaenV1751Board1Label",
                                               "CaenV1751Board1");
  fCaenV1751Board2Label = p.get< std::string >("CaenV1751Board2Label",
                                               "CaenV1751Board2");
}

//------------------------------------------------------------------------------
void DigitReader::beginJob()
{
  return;
}

//------------------------------------------------------------------------------
void DigitReader::analyze(art::Event const & evt)
{

  std::cout << "evt.run(): " << evt.run() << "; evtsubRun(): " << evt.subRun()
            << "; evt.event(): " << evt.event() << std::endl;

  //std::cout << "evt.event(): " << evt.event() << std::endl;

  art::Handle< std::vector<raw::AuxDetDigit> > digitVecHandle;
  evt.getByLabel("FragmentToDigit", fCaenV1751Board1Label, digitVecHandle);
  std::cout << "digitVecHandle Size is: " << digitVecHandle->size() << std::endl;

  for (size_t i = 0; i < digitVecHandle->size(); ++i) {

    art::Ptr<raw::AuxDetDigit> digitVec(digitVecHandle, i);
    std::cout << "digitVec->NADC(): " << digitVec->NADC() << std::endl;
    std::cout << "digitVec->Channel(): " << digitVec->Channel() << std::endl;
    std::cout << "digitVec->AuxDetName(): " << digitVec->AuxDetName() << std::endl;
    std::cout << "digitVec->TimeStamp(): " << digitVec->TimeStamp() << std::endl;

  }

  //art::Handle< std::vector<raw::AuxDetDigit> > digitVecHandleDummy;
  //evt.getByLabel("FragmentToDigit", "a", digitVecHandleDummy);
  //std::cout << "digitVecHandleDummy Size is: " << digitVecHandleDummy->size() << std::endl;

  //for (size_t dummyIter = 0; dummyIter < digitVecHandleDummy->size(); ++dummyIter){ // ++ move

  //  art::Ptr<raw::AuxDetDigit> dummydigitVec(digitVecHandleDummy, dummyIter);
  //  std::cout << "dummydigitVec->NADC(): " << dummydigitVec->NADC() << std::endl;
  //  std::cout << "dummydigitVec->Channel(): " << dummydigitVec->Channel() << std::endl;
  //  std::cout << "dummydigitVec->AuxDetName(): " << dummydigitVec->AuxDetName() << std::endl;
  //}

}

DEFINE_ART_MODULE(DigitReader)
