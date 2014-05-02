// HOW TO ACCESS MYSQL FROM C++

// You must have installed the MySQL - C++ connector
// On debian system it is called libmysqlcppconn-dev
// compile this with "g++ -lmysqlcppconn mysqlconn.C"

#include <iostream>
#include "mysql_connection.h"
#include "mysql_driver.h"
#include "/usr/include/cppconn/resultset.h"
#include "/usr/include/cppconn/statement.h"

using namespace std;
using namespace sql;
using namespace sql::mysql;

int main(){

  MySQL_Driver *driver;
  Connection *con;

  driver = sql::mysql::get_mysql_driver_instance();
  con = driver->connect("tcp://hostname:3306", "user", "password");

  sql::Statement *stmt;
  sql::ResultSet  *res;
  stmt = con->createStatement();
  stmt->execute("USE rundb_v1");

  res = stmt->executeQuery("SELECT * FROM run");

  while (res->next()) {
    cout << "run_number = " << res->getInt("run_number");
    cout << ", run_start_user_comment = '" << res->getString("run_start_user_comment") << "'" << endl;
  }

  delete res;
  delete stmt;
  delete con;

  return 0;

}
