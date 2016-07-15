////////////////////////////////////////////////////////////////////////
// \file DetectorPropertiesLArIAT.h
//
// \brief Experiment-specific service to contain information about detector
//  electronics, etc
//
// \author jpaley@fnal.gov
// 
////////////////////////////////////////////////////////////////////////
#ifndef DETECTOR_PROPERTIES_LARIAT_H
#define DETECTOR_PROPERTIES_LARIAT_H

// LArSoft libraries
#include "larcore/Geometry/GeometryCore.h"
#include "larcore/CoreUtils/ProviderPack.h"
#include "lardata/DetectorInfoServices/LArPropertiesService.h"
#include "lardata/DetectorInfo/DetectorClocks.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"

// framework libraries
#include "fhiclcpp/ParameterSet.h"
#include "fhiclcpp/types/Sequence.h"
#include "fhiclcpp/types/Table.h"
#include "fhiclcpp/types/Atom.h"

// C/C++ standard libraries
#include <set>

///General LArSoft Utilities
namespace ldp{
  
  class DetectorPropertiesLArIAT : public detinfo::DetectorProperties {
    public:
      /// List of service providers we depend on
      using providers_type = lar::ProviderPack<
        geo::GeometryCore,
        detinfo::LArProperties,
        detinfo::DetectorClocks
        >;
       
      /// Structure for configuration parameters
      struct Configuration_t {
        using Name = fhicl::Name;
        using Comment = fhicl::Comment;
        
        fhicl::Sequence<double> Efield { Name("Efield"), Comment(
          "electric field in front of each wire plane (the last one is the big one!) [kV/cm]")
          };
    
        fhicl::Atom<double      > Electronlifetime         {
          Name("Electronlifetime"        ),
          Comment("electron lifetime in liquid argon [us]")
        };
        fhicl::Atom<bool        > GetElectronlifetimeFromDB {
          Name("GetElectronlifetimeFromDB"),
          Comment("option to get electron lifetime from LArIAT conditions database")
        };
        fhicl::Atom<double      > Temperature              {
          Name("Temperature"             ),
          Comment("argon temperature [K]")
        };
        fhicl::Atom<double      > ElectronsToADC           {
          Name("ElectronsToADC"          ),
          Comment("conversion factor: (ADC counts)/(ionization electrons)")
        };
        fhicl::Atom<unsigned int> NumberTimeSamples        {
          Name("NumberTimeSamples"       ),
          Comment("number of TPC readout TDC clock ticks per event")
        };
        fhicl::Atom<unsigned int> ReadOutWindowSize        {
          Name("ReadOutWindowSize"       ),
          Comment("number of TPC readout TDC clock ticks per readout window")
        };
        fhicl::Atom<double      > TimeOffsetU              {
          Name("TimeOffsetU"             ),
          Comment("tick offset subtracted to to convert spacepoint coordinates to hit times on view U")
        };
        fhicl::Atom<double      > TimeOffsetV              {
          Name("TimeOffsetV"             ),
          Comment("tick offset subtracted to to convert spacepoint coordinates to hit times on view V")
        };
        fhicl::Atom<double      > TimeOffsetZ              {
          Name("TimeOffsetZ"             ),
          Comment("tick offset subtracted to to convert spacepoint coordinates to hit times on view Z")
        };
        
        fhicl::Atom<double      > SternheimerA             {
          Name("SternheimerA"),
          Comment("parameter a of Sternheimer correction delta = 2log(10) x - cbar + { a (x1-x)^k } theta(x1-x), x = log10(p/m)")
        };
        fhicl::Atom<double      > SternheimerK             {
          Name("SternheimerK"),
          Comment("parameter k of Sternheimer correction delta = 2log(10) x - cbar + { a (x_1-x)^k } theta(x1-x), x = log10(p/m)")
        };
        fhicl::Atom<double      > SternheimerX0            {
          Name("SternheimerX0"),
          Comment("minimum x = log10(p/m) for the application of Sternheimer correction")
        };
        fhicl::Atom<double      > SternheimerX1            {
          Name("SternheimerX1"),
          Comment("parameter x_1 of Sternheimer correction delta = 2log(10) x - cbar + { a (x_1-x)^k } theta(x1-x), x = log10(p/m)")
        };
        fhicl::Atom<double      > SternheimerCbar             {
          Name("SternheimerCbar"),
          Comment("parameter cbar of Sternheimer correction delta = 2log(10) x - cbar + { a (x_1-x)^k } theta(x1-x), x = log10(p/m)")
        };
      
      }; // Configuration_t
 
      DetectorPropertiesLArIAT();
      DetectorPropertiesLArIAT(fhicl::ParameterSet const& pset, 
			 const geo::GeometryCore* geo,
			 const detinfo::LArProperties* lp,
			 const detinfo::DetectorClocks* c,
			 std::set<std::string> ignore_params = {}
			 );
      /**
       * @brief Constructs the provider and sets up the dependencies
       * @param pset FHiCL parameter set for provider configuration
       * @param providers pack of providers DetectorPropertiesLArIAT depends on
       * @see Setup()
       */
      DetectorPropertiesLArIAT(fhicl::ParameterSet const& pset,
                         providers_type providers,
                         std::set<std::string> ignore_params = {});
      DetectorPropertiesLArIAT(DetectorPropertiesLArIAT const&) = delete;
      virtual ~DetectorPropertiesLArIAT() = default;
      
      /**
       * @brief Configures the provider, first validating the configuration
       * @param p configuration parameter set
       * @param ignore_params parameters to be ignored (optional)
       * 
       * This method will validate the parameter set (except for the parameters
       * it's explicitly told to ignore) and extract the useful information out
       * of it.
       */
      void ValidateAndConfigure
        (fhicl::ParameterSet const& p, std::set<std::string> ignore_params = {});

      
      /// Extracts the relevant configuration from the specified object
      void Configure(Configuration_t const& config);
      
      /**
       * @brief Validates the specified configuration
       * @param p configuration parameter set
       * @param ignore_params parameters to be ignored (optional)
       * @return a parsed configuration object
       * @see ValidateAndConfigure(), Configure()
       * 
       * This method will validate the parameter set (except for the parameters
       * it's explicitly told to ignore) and it returns an object ready to
       * be used with Configure().
       */
      Configuration_t ValidateConfiguration
        (fhicl::ParameterSet const& p, std::set<std::string> ignore_params = {})
        ;

      bool Update(uint64_t ts);
      bool UpdateElectronLifetime(uint64_t ts);
      bool UpdateClocks(const detinfo::DetectorClocks* clks);
      
      /**
       * @brief Sets all the providers at once
       * @param providers the pack of service providers we depend on
       * 
       * Example:
       *     
       *     lar::DetectorPropertiesLArIAT::providers_type providers;
       *     providers.set(lar::providerFrom<geo::Geometry>());
       *     providers.set(lar::providerFrom<detinfo::LArPropertiesService>());
       *     providers.set(lar::providerFrom<detinfo::DetectorClocksService>());
       *     detprop->Setup(providers);
       *
       */
      void Setup(providers_type providers);
        
      void SetGeometry(const geo::GeometryCore* g) { fGeo = g; }
      void SetLArProperties(const detinfo::LArProperties* lp) { fLP = lp; }
      void SetDetectorClocks(const detinfo::DetectorClocks* clks) { fClocks = clks; }

      void SetNumberTimeSamples(unsigned int nsamp) { fNumberTimeSamples=nsamp;}
      // Accessors.

      virtual double Efield(unsigned int planegap=0) const override; ///< kV/cm

      virtual double DriftVelocity(double efield=0., double temperature=0.) const override;  ///< cm/us
      
      /// dQ/dX in electrons/cm, returns dE/dX in MeV/cm.
      virtual double BirksCorrection(double dQdX) const override;
      virtual double ModBoxCorrection(double dQdX) const override;

      virtual double ElectronLifetime() const override { return fElectronlifetime; }   //< microseconds
      
      
      /**
       * @brief Returns argon density at a given temperature
       * @param temperature the temperature in kelvin
       * @return argon density in g/cm^3
       * 
       * Density is nearly a linear function of temperature.
       * See the NIST tables for details
       * Slope is between -6.2 and -6.1, intercept is 1928 kg/m^3.
       * This parameterization will be good to better than 0.5%.
       */
      virtual double Density(double temperature) const override;                          ///< g/cm^3
      
      // need to provide a definition, since the override above hides the inherited one
      virtual double Density() const override { return Density(Temperature()); }
      
      /// In kelvin.
      virtual double Temperature()                   const override { return fTemperature; }
      
      /**
       * @brief Restricted mean energy loss (dE/dx)
       * @param mom  momentum of incident particle [GeV/c]
       * @param mass mass of incident particle [GeV/c^2]
       * @param tcut maximum kinetic energy of delta rays [MeV]; 0 for unlimited
       * @return the restricted mean energy loss (dE/dx) in units of MeV/cm
       *
       * Returned value is always positive.
       * For unrestricted mean energy loss, set tcut = 0 (special case),
       * or tcut large.
       * 
       * Based on Bethe-Bloch formula as contained in particle data book.
       * Material parameters are from the configuration.
       */
      virtual double Eloss(double mom, double mass, double tcut) const override;
      
      /**
       * @brief Energy loss fluctuation (@f$ \sigma_{E}^2 / x @f$)
       * @param mom  momentum of incident particle in [GeV/c]
       * @return energy loss fluctuation in MeV^2/cm
       *
       * Based on Bichsel formula referred to but not given in pdg.
       */
      virtual double ElossVar(double mom, double mass) const override;

      virtual double       SamplingRate()      const override { return fTPCClock.TickPeriod() * 1.e3; }
      virtual double       ElectronsToADC()    const override { return fElectronsToADC; }
      virtual unsigned int NumberTimeSamples() const override { return fNumberTimeSamples; }
      virtual unsigned int ReadOutWindowSize() const override { return fReadOutWindowSize; }
      virtual int          TriggerOffset()     const override;
      virtual double       TimeOffsetU()       const override{ return fTimeOffsetU; };
      virtual double       TimeOffsetV()       const override { return fTimeOffsetV; };
      virtual double       TimeOffsetZ()       const override{ return fTimeOffsetZ; };

      virtual double       ConvertXToTicks(double X, int p, int t, int c) const override;
      virtual double       ConvertXToTicks(double X, geo::PlaneID const& planeid) const override
        { return ConvertXToTicks(X, planeid.Plane, planeid.TPC, planeid.Cryostat); }
      virtual double       ConvertTicksToX(double ticks, int p, int t, int c) const override;
      virtual double       ConvertTicksToX(double ticks, geo::PlaneID const& planeid) const override
        { return ConvertTicksToX(ticks, planeid.Plane, planeid.TPC, planeid.Cryostat); }

      virtual double       GetXTicksOffset(int p, int t, int c) const override;
      virtual double       GetXTicksOffset(geo::PlaneID const& planeid) const override
        { return GetXTicksOffset(planeid.Plane, planeid.TPC, planeid.Cryostat); }
      virtual double       GetXTicksCoefficient(int t, int c) const override;
      virtual double       GetXTicksCoefficient(geo::TPCID const& tpcid) const override
        { return GetXTicksCoefficient(tpcid.TPC, tpcid.Cryostat); }
      virtual double       GetXTicksCoefficient() const override;

      // The following methods convert between TDC counts (SimChannel time) and
      // ticks (RawDigit/Wire time).
      virtual double       ConvertTDCToTicks(double tdc) const override;
      virtual double       ConvertTicksToTDC(double ticks) const override;
      
      
      /// Verifies that the provider is in a fully configured status
      /// @throw cet::exception (category DetectorPropertiesLArIAT) if not ok
      void CheckIfConfigured() const;
      
    protected:
      
     
      /// Parameters for Sternheimer density effect corrections
      struct SternheimerParameters_t {
        double a;               ///< parameter a
        double k;               ///< parameter k
        double x0;              ///< parameter x0
        double x1;              ///< parameter x1
        double cbar;            ///< parameter Cbar
      }; //  SternheimerParameters_t
      
      void         CalculateXTicksParams();
      
      // service providers we depend on;
      // in principle could be replaced by a single providerpacl_type.
      const detinfo::LArProperties* fLP;
      const detinfo::DetectorClocks* fClocks;
      const geo::GeometryCore* fGeo;

      bool  fGetElectronlifetimeFromDB;
      
      std::vector< double >          fEfield;           ///< kV/cm (per inter-plane volume)
      double                         fElectronlifetime; ///< microseconds
      double                         fTemperature;      ///< kelvin
      double       fSamplingRate;      ///< in ns
      double 	   fElectronsToADC;    ///< conversion factor for # of ionization electrons to 1 ADC count
      unsigned int fNumberTimeSamples; ///< number of clock ticks per event
      unsigned int fReadOutWindowSize; ///< number of clock ticks per readout window
      double       fTimeOffsetU;       ///< time offset to convert spacepoint coordinates to hit times on view U
      double       fTimeOffsetV;       ///< time offset to convert spacepoint coordinates to hit times on view V
      double       fTimeOffsetZ;       ///< time offset to convert spacepoint coordinates to hit times on view Z
      
      SternheimerParameters_t fSternheimerParameters; ///< Sternheimer parameters
      
      double       fXTicksCoefficient; ///< Parameters for x<-->ticks

      std::vector<std::vector<std::vector<double> > > fXTicksOffsets;
      std::vector<std::vector<double> >               fDriftDirection;

      ::detinfo::ElecClock fTPCClock;     ///< TPC electronics clock
    }; // class DetectorPropertiesLArIAT
} //namespace detinfo

#endif // DETINFO_DETECTOR_PROPERTIES_H
