PS C:\Users\Tim\Documents\ETH Sem5\ETH Sem5_CAstro\nBodyProject> g++ tests2.cpp -o tests2.exe -I C:/Python38/include -I C:/Python38/Lib/site-packages/numpy/core/include -L C:/Python38/libs -lpython38 
In file included from C:/Python38/include/Python.h:63,
                 from matplotlibcpp.h:5,
                 from tests2.cpp:7:
C:/Python38/include/pyport.h:717: warning: "LONG_BIT" redefined
  717 | #define LONG_BIT (8 * SIZEOF_LONG)
      |
In file included from /usr/lib/gcc/x86_64-pc-cygwin/11/include/limits.h:203,
                 from /usr/lib/gcc/x86_64-pc-cygwin/11/include/syslimits.h:7,
                 from /usr/lib/gcc/x86_64-pc-cygwin/11/include/limits.h:34,
                 from C:/Python38/include/Python.h:11,
                 from matplotlibcpp.h:5,
                 from tests2.cpp:7:
/usr/include/limits.h:31: note: this is the location of the previous definition
   31 | #define LONG_BIT (__SIZEOF_LONG__ * __CHAR_BIT__)
      |
In file included from /usr/include/sys/stat.h:22,
                 from C:/Python38/include/pyport.h:245,
                 from C:/Python38/include/Python.h:63,
                 from matplotlibcpp.h:5,
                 from tests2.cpp:7:
C:/Python38/include/fileutils.h:85:12: error: expected ‘;’ at end of member declaration
   85 |     time_t st_atime;
      |            ^~~~~~~~
C:/Python38/include/fileutils.h:85:12: error: expected unqualified-id before ‘.’ token
   85 |     time_t st_atime;
      |            ^~~~~~~~
C:/Python38/include/fileutils.h:87:12: error: expected ‘;’ at end of member declaration
   87 |     time_t st_mtime;
      |            ^~~~~~~~
C:/Python38/include/fileutils.h:87:12: error: expected unqualified-id before ‘.’ token
   87 |     time_t st_mtime;
      |            ^~~~~~~~
C:/Python38/include/fileutils.h:89:12: error: expected ‘;’ at end of member declaration
   89 |     time_t st_ctime;
      |            ^~~~~~~~
C:/Python38/include/fileutils.h:89:12: error: expected unqualified-id before ‘.’ token
   89 |     time_t st_ctime;