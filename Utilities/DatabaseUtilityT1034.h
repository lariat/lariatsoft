//////////////////////////////////////////////////////////////
// Name:      DatabaseUtilityT1034.h
// Date:      10 July 2015
// Author:    Everybody is an author!
//////////////////////////////////////////////////////////////
// DatabaseUtilityT1034 class for LArIAT.
// snake_case is better than CamelCase.
//////////////////////////////////////////////////////////////

#ifndef DATABASEUTILITYT1034_H
#define DATABASEUTILITYT1034_H

// Framework includes
#include "art/Framework/Services/Registry/ActivityRegistry.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Registry/ServiceMacros.h"
#include "fhiclcpp/ParameterSet.h"

// PostgreSQL includes
#include <libpq-fe.h>

namespace util {

  //-----------------------------------------------------------------------
  //-----------------------------------------------------------------------
  // class definition

  class DatabaseUtilityT1034 {

   public:

    // standard constructor and destructor for an art service
    DatabaseUtilityT1034(fhicl::ParameterSet   const& pset,
                         art::ActivityRegistry      & reg);
    ~DatabaseUtilityT1034();

    // this method reads in any parameters from the .fcl files
    void reconfigure(fhicl::ParameterSet const& pset);

    // retrieve a single value with a single query
    std::string GetValue(std::string const& Query);

    // retrieve a row of values with a single query
    // returns a map where the key is the column name and the value is the value
    std::map< std::string, std::string > GetValues(std::string const& Query);

    // retrieve a single value from the lariat_xml_database table
    std::string GetConfigValue(std::string const& ColumnName,
                               int         const& RunNumber);

    // retrieve a single value from the lariat_ifbeam_database table
    std::string GetIFBeamValue(std::string const& ColumnName,
                               int         const& RunNumber,
                               int         const& SubRunNumber);

    //////////////////////////////////////////////////////////
    // the following methods return a map where the key is
    // the column name and the value is the value
    //////////////////////////////////////////////////////////

    // retrieve row from the lariat_xml_database table
    std::map< std::string, std::string > GetConfigValues(std::vector< std::string > const& ColumnName,
                                                         int                        const& RunNumber);

    // retrieve row from the lariat_ifbeam_database table
    std::map< std::string, std::string > GetIFBeamValues(std::vector< std::string > const& ColumnName,
                                                         int                        const& RunNumber,
                                                         int                        const& SubRunNumber);

    // retrieve entire row from the lariat_xml_database table
    std::map< std::string, std::string > GetAllConfigValues(int const& RunNumber);

    // retrieve entire row from the lariat_ifbeam_database table
    std::map< std::string, std::string > GetAllIFBeamValues(int const& RunNumber,
                                                            int const& SubRunNumber);

    //////////////////////////////////////////////////////////

    // this method is used for testing porpoises
    void hello_world();

   private:

    // this method is used for testing dolphins
    void hello_kitty();

    bool Connect();
    bool Disconnect();

    std::string fDBHost;
    std::string fDBPort;
    std::string fDBName;
    std::string fDBUser;
    std::string fDBPasswordFile;
    std::string fDBPassword;
    unsigned int fDBReconnectWaitTime;
    unsigned int fDBNumberConnectAttempts;

    std::string fConfigTableName;
    std::string fIFBeamTableName;
    std::string fLHCdbTableName;

    std::string fConnectionInfo;

    PGconn * fConnection;

  };

}

DECLARE_ART_SERVICE(util::DatabaseUtilityT1034, LEGACY)

#endif
