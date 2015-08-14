//////////////////////////////////////////////////////////////
// Name:      DatabaseUtilityT1034_service.cc
// Date:      10 July 2015
// Author:    Everybody is an author!
//////////////////////////////////////////////////////////////
// DatabaseUtilityT1034 service for LArIAT.
//////////////////////////////////////////////////////////////

// Class include
#include "Utilities/DatabaseUtilityT1034.h" 

// Framework includes
#include "messagefacility/MessageLogger/MessageLogger.h"

// C++ includes
#include <iostream>
#include <fstream>

namespace util {

  //-----------------------------------------------------------------------
  //-----------------------------------------------------------------------
  // class implementation

  //-----------------------------------------------------------------------
  // constructor
  DatabaseUtilityT1034::DatabaseUtilityT1034(fhicl::ParameterSet   const& pset,
                                             art::ActivityRegistry      & reg) {

    // read in parameters from .fcl files
    this->reconfigure(pset);

  }

  //-----------------------------------------------------------------------
  // destructor
  DatabaseUtilityT1034::~DatabaseUtilityT1034() {}

  //-----------------------------------------------------------------------
  void DatabaseUtilityT1034::reconfigure(fhicl::ParameterSet const& pset) {

    // get parameters from .fcl files
    fDBHost = pset.get< std::string >("DBHost", "ifdb02.fnal.gov");
    fDBPort = pset.get< std::string >("DBPort", "5443");
    fDBName = pset.get< std::string >("DBName", "lariat_prd");
    fDBUser = pset.get< std::string >("DBUser", "lariat_prd_user");

    fDBPasswordFile = pset.get< std::string >("DBPasswordFile", "Utilities/lariat_prd_passwd");
    fDBReconnectWaitTime = pset.get< unsigned int >("DBReconnectWaitTime", 10);
    fDBNumberConnectAttempts = pset.get< unsigned int >("DBNumberConnectAttempts", 3);

    fConfigTableName = pset.get< std::string >("ConfigTableName", "lariat_xml_database");
    fIFBeamTableName = pset.get< std::string >("IFBeamTableName", "lariat_ifbeam_database");

    // find the password file
    cet::search_path search_path("FW_SEARCH_PATH");
    std::string password_file_path;
    search_path.find_file(fDBPasswordFile, password_file_path);

    // read in the password from the password file
    std::ifstream in(password_file_path.c_str());
    if (in.is_open()) {
      std::getline(in, fDBPassword);
      in.close();
    }

    // set connection info for connecting to database
    fConnectionInfo = "host = "      + fDBHost +
                      " port = "     + fDBPort +
                      " dbname = "   + fDBName +
                      " user = "     + fDBUser +
                      " password = " + fDBPassword;

  }

  //-----------------------------------------------------------------------
  void DatabaseUtilityT1034::hello_world() {

    mf::LogVerbatim("DatabaseUtilityT1034")
        << "\n///////////////////////////////////////////////////////////\n"
        << "Hello, World!\n"
        << "///////////////////////////////////////////////////////////\n";

    this->hello_kitty();

    return;
  }

  //-----------------------------------------------------------------------
  void DatabaseUtilityT1034::hello_kitty() {

    bool connected = this->Connect();

    if (!connected) return;

    std::string query = "select count(*) from public.lariat_xml_database;";

    PGresult * result;
    result = PQexec(fConnection, query.c_str());

    int number_fields = PQnfields(result);
    int number_rows = PQntuples(result);

    mf::LogVerbatim("DatabaseUtilityT1034") << "Query: " << query;

    for (int i = 0; i < number_rows; ++i) {
      for (int j = 0; j < number_fields; ++j) {
        mf::LogVerbatim("DatabaseUtilityT1034")
            << "Column: " << PQfname(result, j) << "; "
            << "Value: " << PQgetvalue(result, i, j);
      }
    }

    PQclear(result);

    this->Disconnect();

    return;
  }

  //-----------------------------------------------------------------------
  bool DatabaseUtilityT1034::Connect() {

    // number of attempts made to connect to database
    size_t number_attempts = 0;
    size_t max_number_attempts = fDBNumberConnectAttempts;

    mf::LogVerbatim("DatabaseUtilityT1034") << "Attempting to connect to database...";

    // create a connection to the database
    fConnection = PQconnectdb(fConnectionInfo.c_str());

    // attempt to reconnect if connection is not OK
    while (PQstatus(fConnection) != CONNECTION_OK) {

      mf::LogVerbatim("DatabaseUtilityT1034") << "Connection to database failed: "
                                              << PQerrorMessage(fConnection);
      PQfinish(fConnection);

      number_attempts += 1;

      // abort if too many attempts are made
      if (number_attempts > max_number_attempts) {

        mf::LogError("DatabaseUtilityT1034")
            << "Too many attempts made to connect to the database... Aborting.";

        return false;
      }

      // wait before reconnecting
      sleep(fDBReconnectWaitTime);

      mf::LogVerbatim("DatabaseUtilityT1034")
          << "Attempting to reconnect to database... attempt number "
          << number_attempts;

      // create a connection to the database
      fConnection = PQconnectdb(fConnectionInfo.c_str());

    }

    mf::LogVerbatim("DatabaseUtilityT1034") << "Connection established!";

    return true;
  }

  //-----------------------------------------------------------------------
  bool DatabaseUtilityT1034::Disconnect() {

    // if there is a connection, close it
    if (fConnection) {
      PQfinish(fConnection);
      mf::LogVerbatim("DatabaseUtilityT1034") << "Connection closed.";
    }

    return true;
  }

  //-----------------------------------------------------------------------
  std::string DatabaseUtilityT1034::GetValue(std::string const& Query) {

    // connect to database
    bool connected = this->Connect();

    // if no connection can be established, return empty string
    if (!connected) return "";

    mf::LogVerbatim("DatabaseUtilityT1034") << "Query: " << Query;

    // query the database, and get the result of query
    PGresult * result;
    result = PQexec(fConnection, Query.c_str());

    // get number of columns and rows of the result
    int number_fields = PQnfields(result);
    int number_rows = PQntuples(result);

    // if the result has more than one column or row, return empty string
    if (number_fields != 1 or number_rows != 1) return "";

    // cast the value of result into a string
    std::string value = std::string(PQgetvalue(result, 0, 0));

    // destroy the evidence
    PQclear(result);

    // fall off the grid and disappear
    this->Disconnect();

    return value;
  }

  //-----------------------------------------------------------------------
  std::map< std::string, std::string > DatabaseUtilityT1034::GetValues(std::string const& Query) {

    std::map< std::string, std::string > values;
    values[""] = "";

    // connect to database
    bool connected = this->Connect();

    // if no connection can be established, return empty string
    if (!connected) return values;

    mf::LogVerbatim("DatabaseUtilityT1034") << "Query: " << Query;

    // query the database, and get the result of query
    PGresult * result;
    result = PQexec(fConnection, Query.c_str());

    // get number of columns and rows of the result
    int number_fields = PQnfields(result);
    int number_rows = PQntuples(result);

    // if the result has more than one row, return empty string
    if (number_rows != 1) return values;

    // fill map with column name as the key and value as the mapped value
    for (int i = 0; i < number_rows; ++i) {
      for (int j = 0; j < number_fields; ++j) {
        //mf::LogVerbatim("DatabaseUtilityT1034")
        //    << "Column: " << PQfname(result, j) << "; "
        //    << "Value: " << PQgetvalue(result, i, j);
        values[std::string(PQfname(result, j))] = std::string(PQgetvalue(result, i, j));
      }
    }

    // destroy the evidence
    PQclear(result);

    // fall off the grid and disappear
    this->Disconnect();

    return values;
  }

  //-----------------------------------------------------------------------
  std::string DatabaseUtilityT1034::GetConfigValue(std::string const& ColumnName,
                                                   int         const& RunNumber) {

    // set schema and table names
    std::string schema_name = "public";
    std::string table_name = fConfigTableName;

    // set query
    std::string query = "select " + ColumnName +
                        " from " + schema_name + "." + table_name +
                        " where runnumber = " + std::to_string(RunNumber) + ";";

    // get value with query
    std::string value = this->GetValue(query);

    return value;
  }

  //-----------------------------------------------------------------------
  std::string DatabaseUtilityT1034::GetIFBeamValue(std::string const& ColumnName,
                                                   int         const& RunNumber,
                                                   int         const& SubRunNumber) {

    // set schema and table names
    std::string schema_name = "public";
    std::string table_name = fIFBeamTableName;

    // set query
    std::string query = "select " + ColumnName +
                        " from " + schema_name + "." + table_name +
                        " where runnumber = " + std::to_string(RunNumber) +
                        " and subrun = " + std::to_string(SubRunNumber) + ";";

    // get value with query
    std::string value = this->GetValue(query);

    return value;
  }

  //-----------------------------------------------------------------------
  std::map< std::string, std::string > DatabaseUtilityT1034::GetConfigValues(std::vector< std::string > const& ColumnNames,
                                                                             int                        const& RunNumber) {

    // concatenate column names into single string for query
    std::string column_names;

    size_t number_columns = ColumnNames.size();

    for (size_t i = 0; i < number_columns; ++i) {
      if (i == 0) {
        column_names = ColumnNames.at(i);
      }
      else {
        column_names += ", " + ColumnNames.at(i);
      }
    }

    // set schema and table names
    std::string schema_name = "public";
    std::string table_name = fConfigTableName;

    // set query
    std::string query = "select " + column_names +
                        " from " + schema_name + "." + table_name +
                        " where runnumber = " + std::to_string(RunNumber) + ";";

    // get values with query
    std::map< std::string, std::string > values = this->GetValues(query);

    return values;
  }

  //-----------------------------------------------------------------------
  std::map< std::string, std::string > DatabaseUtilityT1034::GetIFBeamValues(std::vector< std::string > const& ColumnNames,
                                                                             int                        const& RunNumber,
                                                                             int                        const& SubRunNumber) {

    // concatenate column names into single string for query
    std::string column_names;

    size_t number_columns = ColumnNames.size();

    for (size_t i = 0; i < number_columns; ++i) {
      if (i == 0) {
        column_names = ColumnNames.at(i);
      }
      else {
        column_names += ", " + ColumnNames.at(i);
      }
    }

    // set schema and table names
    std::string schema_name = "public";
    std::string table_name = fIFBeamTableName;

    // set query
    std::string query = "select " + column_names +
                        " from " + schema_name + "." + table_name +
                        " where runnumber = " + std::to_string(RunNumber) +
                        " and subrun = " + std::to_string(SubRunNumber) + ";";

    // get values with query
    std::map< std::string, std::string > values = this->GetValues(query);

    return values;
  }

  //-----------------------------------------------------------------------
  std::map< std::string, std::string > DatabaseUtilityT1034::GetAllConfigValues(int const& RunNumber) {

    // set schema and table names
    std::string schema_name = "public";
    std::string table_name = fConfigTableName;

    // set query
    std::string query = "select * from " + schema_name + "." + table_name +
                        " where runnumber = " + std::to_string(RunNumber) + ";";

    // get values with query
    std::map< std::string, std::string > values = this->GetValues(query);

    return values;
  }

  //-----------------------------------------------------------------------
  std::map< std::string, std::string > DatabaseUtilityT1034::GetAllIFBeamValues(int const& RunNumber,
                                                                                int const& SubRunNumber) {

    // set schema and table names
    std::string schema_name = "public";
    std::string table_name = fIFBeamTableName;

    // set query
    std::string query = "select * from " + schema_name + "." + table_name +
                        " where runnumber = " + std::to_string(RunNumber) +
                        " and subrun = " + std::to_string(SubRunNumber) + ";";

    // get values with query
    std::map< std::string, std::string > values = this->GetValues(query);

    return values;
  }

}

namespace util {

  DEFINE_ART_SERVICE(DatabaseUtilityT1034)

} // namespace util
