/**
 * @file   AuxDetGeometryCore.h
 * @brief  Access the description of auxiliary detector geometry
 * @author brebel@fnal.gov
 * @see    AuxDetGeometryCore.cxx
 *
 * Structure of the header:
 *     
 *     namespace geo {
 *       
 *       // forward class declarations
 *       
 *       namespace details {
 *         
 *         // geometry iterator base class
 *       
 *       }
 *       
 *       // geometry iterators declaration
 *       //  - cryostat_id_iterator
 *       //  - TPC_id_iterator
 *       //  - plane_id_iterator
 *       //  - wire_id_iterator
 *       
 *       // GeometryData_t definition (part of GeometryCore)
 *       
 *       // GeometryCore declaration
 *     
 *     }
 *     
 *
 *
 * Revised <seligman@nevis.columbia.edu> 29-Jan-2009
 *         Revise the class to make it into more of a general detector interface
 * Revised <petrillo@fnal.gov> 27-Apr-2015
 *         Factorization into a framework-independent GeometryCore.h and a
 *         art framework interface
 * Revised <petrillo@fnal.gov> 30-Apr-2015
 *         Redesign of the iterators
 */
#ifndef GEO_AUXDETGEOMETRYCORE_H
#define GEO_AUXDETGEOMETRYCORE_H

// LArSoft libraries

// Framework and infrastructure libraries
#include "fhiclcpp/ParameterSet.h"

// ROOT libraries
#include <TVector3.h>

// C/C++ standard libraries
#include <cstddef> // size_t
#include <string>
#include <vector>
#include <set>
#include <memory> // std::shared_ptr<>
#include <iterator> // std::forward_iterator_tag
#include <type_traits> // std::is_base_of<>

#include "Geo/AuxDetChannelMapAlg.h"

// ROOT class prototypes
class TGeoManager;
class TGeoNode;
class TGeoMaterial;


/// Namespace collecting geometry-related classes utilities
namespace geo {
  
  
  // Forward declarations within namespace.
  class AuxDetGeo;
  class AuxDetSensitiveGeo;
  class AuxDetGeometryCore;
  
  
  /// Data in the geometry description
  struct AuxDetGeometryData_t {
    
    /// Type of list of auxiliary detectors
    using AuxDetList_t = std::vector<AuxDetGeo*>;
    
    AuxDetList_t   auxDets;   ///< The auxiliary detectors
    
  }; // AuxDetGeometryData_t
  
  
  /** **************************************************************************
   * @brief Description of geometry of one set of auxiliary detectors
   * 
   * @note All lengths are specified in centimetres
   * 
   * 
   * How to correctly instantiate a GeometryCore object
   * ---------------------------------------------------
   * 
   * Instantiation is a multi-step procedure:
   * 1. construct a GeometryCore object (the "service provider"),
   *    with the full configuration; at this step, configuration is just stored
   * 2. load a geometry with GeometryCore::LoadGeometryFile();
   *    this loads the detector geometry information
   * 3. prepare a channel map algorithm object (might use for example
   *    GeometryCore::DetectorName() or the detector geometry from the
   *    newly created object, but any use of channel mapping related functions
   *    is forbidden and it would yield undefined behaviour (expected to be
   *    catastrophic)
   * 4. acquire the channel mapping algorithm with
   *    GeometryCore::ApplyChannelMap(); at this point, the ChannelMapAlg object
   *    is asked to initialize itself and to perform whatever modifications to
   *    the geometry provider is needed.
   * 
   * Step 3 (creation of the channel mapping algorithm object) can be performed
   * at any time before step 4, provided that no GeometryCore instance is needed
   * for it.
   * 
   * 
   * Configuration parameters
   * -------------------------
   * 
   * - *Name* (string; mandatory): string identifying the detector; it can be
   *   different from the base name of the file used to initialize the geometry;
   *   standard names are recommended by each experiment.
   *   This name can be used, for example, to select which channel mapping
   *   algorithm to use.
   * - *SurfaceY* (real; mandatory): depth of the detector, in centimetrs;
   *   see SurfaceY() for details
   * - *MinWireZDist* (real; default: 3)
   * - *PositionEpsilon* (real; default: 0.01%) set the default tolerance
   *   (see DefaultWiggle())
   * 
   */
  class AuxDetGeometryCore {
  public:
    
    /// Type of list of auxiliary detectors
    using AuxDetList_t = AuxDetGeometryData_t::AuxDetList_t;    
    
    /**
     * @brief Initialize geometry from a given configuration
     * @param pset configuration parameters
     * 
     * This constructor does not load any geometry description.
     * The next step is to do exactly that, by GeometryCore::LoadGeometryFile().
     */
    AuxDetGeometryCore(fhicl::ParameterSet const& pset);
    
    /// Destructor
    ~AuxDetGeometryCore();
   
    // You shall not copy or move or assign me!
    AuxDetGeometryCore(AuxDetGeometryCore const&) = delete;
    AuxDetGeometryCore(AuxDetGeometryCore&&) = delete;
    AuxDetGeometryCore& operator= (AuxDetGeometryCore const&) = delete;
    AuxDetGeometryCore& operator= (AuxDetGeometryCore&&) = delete;
    
    
    /**
     * @brief Returns the full directory path to the geometry file source
     * @return the full directory path to the geometry file source
     * 
     * This is the full path of the source of the detector geometry GeometryCore
     * relies on.
     */
    std::string ROOTFile() const { return fROOTfile; }

    /**
     * @brief Returns the full directory path to the GDML file source
     * @return the full directory path to the GDML file source
     * 
     * This is the full path of the source of the detector geometry handed to
     * the detector simulation (GEANT).
     */
    std::string GDMLFile() const { return fGDMLfile; }
            

    /// Returns a string with the name of the detector, as configured
    std::string DetectorName() const { return fDetectorName; }    
    
    //
    // object description and information
    //
        
    /// @todo use a AutDetID_t instead of unsigned int?
    
    //
    // group features
    //
    
    /**
     * @brief Returns the number of auxiliary detectors
     * 
     * This method returns the total number of scintillator paddles
     * (Auxiliary Detectors aka AuxDet) outside of the cryostat
     * 
     * @todo Change return type to size_t
     */
    unsigned int NAuxDets() const { return AuxDets().size(); }
    
    /**
     * @brief Returns the number of sensitive components of auxiliary detector
     * @param aid ID of the auxiliary detector
     * @return number of sensitive components in the auxiliary detector aid
     * @thrws cet::exception (category "Geometry") if aid does not exist
     */
    unsigned int NAuxDetSensitive(size_t const& aid) const;
    
    //
    // access
    //
    
    /// Returns the full list of pointer to the auxiliary detectors
    std::vector<AuxDetGeo*> const& AuxDetGeoVec() const { return AuxDets(); }
    
    /**
     * @brief Returns the specified auxiliary detector
     * @param ad the auxiliary detector index
     * @return a constant reference to the specified auxiliary detector
     * 
     * @todo what happens if it does not exist?
     * @todo remove the default parameter?
     */
    AuxDetGeo const& AuxDet(unsigned int const ad = 0) const;
    
    /**
     * @brief Returns the index of the auxiliary detector at specified location
     * @param worldLoc 3D coordinates of the point (world reference frame)
     * @return the index of the detector, or UINT_MAX if no detector is there
     * 
     * @todo replace with numeric_limits<>?
     */
    unsigned int FindAuxDetAtPosition(double const worldLoc[3]) const;
    
    /**
     * @brief Fills the indices of the sensitive auxiliary detector at location
     * @param worldLoc 3D coordinates of the point (world reference frame)
     * @param adg (output) auxiliary detector index
     * @param sv (output) sensitive volume index
     */
    void  FindAuxDetSensitiveAtPosition(double const worldLoc[3],
                                        size_t     & adg,
                                        size_t     & sv) const;
    
    /**
     * @brief Returns the auxiliary detector at specified location
     * @param worldLoc 3D coordinates of the point (world reference frame)
     * @param ad (output) the auxiliary detector index
     * @return constant reference to AuxDetGeo object of the auxiliary detector
     * 
     * @todo what happens if it does not exist?
     */
    AuxDetGeo const& PositionToAuxDet(double const worldLoc[3], 
				      unsigned int &ad) const;
    
    /**
     * @brief Returns the auxiliary detector at specified location
     * @param worldLoc 3D coordinates of the point (world reference frame)
     * @param ad (output) the auxiliary detector index
     * @param sv (output) the auxiliary detector sensitive volume index
     * @return reference to AuxDetSensitiveGeo object of the auxiliary detector
     * 
     * @todo what happens if it does not exist?
     */
    const AuxDetSensitiveGeo& PositionToAuxDetSensitive(double const worldLoc[3],
                                                        size_t     & ad,
                                                        size_t     & sv) const;
    
    const uint32_t           PositionToAuxDetChannel(double const worldLoc[3],
						     size_t     & ad,
						     size_t     & sv) const;
    const TVector3           AuxDetChannelToPosition(uint32_t    const& channel,
						     std::string const& auxDetName) const;


    const AuxDetGeo&         ChannelToAuxDet(std::string const& auxDetName,
					     uint32_t    const& channel) const; // return the AuxDetGeo for the given detector 
                                                                                // name and channel

    const AuxDetSensitiveGeo& ChannelToAuxDetSensitive(std::string const& auxDetName,
						       uint32_t    const& channel) const; // return the AuxDetSensitiveGeo for the given
    
    /// @name Geometry initialization
    /// @{
    
    /**
     * @brief Loads the geometry information from the specified files
     * @param gdmlfile path to file to be used for Geant4 simulation
     * @param rootfile path to file for internal geometry representation
     * @see ApplyChannelMap()
     *
     * Both paths must directly resolve to an available file, as no search
     * is performed for them.
     * 
     * The gdmlfile parameter does not have to necessarily be in GDML format,
     * as long as it's something supported by Geant4. This file is not used by
     * the geometry, but its path is provided on request by the simulation
     * modules (see LArSoft `LArG4` module).
     * The rootfile also does not need to be a ROOT file, but just anything
     * that TGeoManager::Import() supports. This file is parsed immediately
     * and the internal geometry representation is built out of it.
     * 
     * @note After calling this method, the detector geometry information can
     * be considered complete, but the geometry service provider is not fully
     * initialized yet, since it's still necessary to provide or update the
     * channel mapping.
     */
    void LoadGeometryFile(std::string gdmlfile, std::string rootfile);
    
    /**
     * @brief Initializes the geometry to work with this channel map
     * @param pChannelMap a pointer to the channel mapping algorithm to be used
     * @see LoadGeometryFile()
     * 
     * The specified channel mapping is used with this geometry.
     * The algorithm object is asked and allowed to make the necessary
     * modifications to the geometry description.
     * These modifications typically involve some resorting of the objects.
     * 
     * The ownership of the algorithm object is shared, usually with a calling
     * framework: we maintain it alive as long as we need it (and no other code
     * can delete it), and we delete it only if no other code is sharing the
     * ownership.
     * 
     * This method needs to be called after LoadGeometryFile() to complete the
     * geometry initialization.
     */
    void ApplyChannelMap(std::shared_ptr<geo::AuxDetChannelMapAlg> pChannelMap);
    /// @}
    
    
  protected:

    /// Returns the object handling the channel map
    geo::AuxDetChannelMapAlg const* AuxDetChannelMap() const { return fChannelMapAlg.get(); }
    
    //@{
    /// Return the internal auxiliary detectors list
    AuxDetList_t&       AuxDets()       { return fGeoData.auxDets; }
    AuxDetList_t const& AuxDets() const { return fGeoData.auxDets; }
    //@}
    
  private:
    
    void FindAuxDet(std::vector<const TGeoNode*>& path, unsigned int depth);
    
    void MakeAuxDet(std::vector<const TGeoNode*>& path, int depth);
    
    /// Deletes the detector geometry structures
    void ClearGeometry();
    
    AuxDetGeometryData_t fGeoData;  ///< The detector description data
    
    std::string    fDetectorName;   ///< Name of the detector.
    std::string    fGDMLfile;       ///< path to geometry file used for Geant4 simulation
    std::string    fROOTfile;       ///< path to geometry file for geometry in GeometryCore
    std::shared_ptr<const geo::AuxDetChannelMapAlg> fChannelMapAlg;  ///< Object containing the channel to wire mapping
  }; // class GeometryCore
  
} // namespace geo




#endif // GEO_AUXDETGEOMETRYCORE_H
