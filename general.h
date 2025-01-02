#ifndef GENERAL_H
#define GENERAL_H
/*============================================================================*/
/*!\file
 * \brief general.h - definitions of general use to (standard) C programs
 *
 * \author  David Hoadley <vcrumble@westnet.com.au>
 *
 * \details
 *      Definitions that almost all of our C programs can find useful. Defines
 *      bool and fixed size int data types, if the compiler standard does not
 *      (i.e. if this is a pre-C99 compiler). Defines a static_assert() macro
 *      if this is a pre-C11 compiler. Some other useful constants and macros
 *      are defined, including macros for assertions in an embedded context.
 * 
 *==============================================================================
 */
/* Which C standard(s) does this compiler support here? */
#if defined(__STDC__)
#  define PREDEF_STANDARD_C_1989
#  if defined(__STDC_VERSION__)
#    define PREDEF_STANDARD_C_1990
#    if (__STDC_VERSION__ >= 199409L)
#      define PREDEF_STANDARD_C_1994
#    endif
#    if (__STDC_VERSION__ >= 199901L)
#      define PREDEF_STANDARD_C_1999
#    endif
#    if (__STDC_VERSION__ >= 201112L)
#      define PREDEF_STANDARD_C_2011
#    endif
#    if (__STDC_VERSION__ >= 201710L)
#      define PREDEF_STANDARD_C_2018
#    endif
#  endif
#endif

#include <cstdint>
#include <cassert>

/*------------------------------------------------------------------------------
 *
 * Peculiar design decisions in the C language
 *
 *----------------------------------------------------------------------------*/

/*! C has some very silly default behaviours. One of these is to define all 
 *  symbols declared at file level as global in scope (i.e. visible to the
 *  linker). This is a recipe for trouble. Secondly, C overloads the 'static'
 *  storage class (which defines how data is allocated to memory) with a second
 *  completely independent meaning - that the scope of the symbol name will be
 *  restricted to the current C source file only. While this is behaviour that
 *  we want, the use of the word 'static' is odd, to say the least. So, use the
 *  following definitions to discipline yourself to sensible behaviour - define
 *  truly global variables and functions as 'GLOBAL', and those which are to be
 *  visible only in the current C source file as 'LOCAL'.
 *  [Each global function definition should be matched with a corresponding
 *   function prototype in the corresponding header file, and each global
 *   variable definition should be matched with a corresponding 'extern' 
 *   declaration in that same header file. There should be no 'extern'
 *   declarations outside of header files.] */  
#define LOCAL   static
#define GLOBAL                  /*!< See above */

/* Another odd design decision was the requirement to insert the keyword "break"
 * after every element of a switch-case construct if you did not want to run on
 * and execute the code of the next element. That is, most of the time you must
 * insert the word "break", and only occasionally do you write code in which you
 * don't want it. So on average, if a "break" is missing in a case statement, it
 * is a bug. So, on those rare occasions when we do intend to run on, we want to
 * let people know that, this time, it isn't a bug. */
#define NOBREAK     /*!< Show reader that 'break' was not omitted by accident */
#define FALLTHROUGH /*!< or if you prefer, use this word instead. */


/*------------------------------------------------------------------------------
 *
 * Useful arithmetic constants
 *
 *----------------------------------------------------------------------------*/

/* A few definitions below are followed by an empty comment containing "!< . "
   and nothing else. This is just to get Doxygen to include the definitions in
   its lists and cross-references. But you don't need me to explain what each
   one is, do you? */
/*      Double precision constants */
#define PI          3.1415926535897932384626433832795028841971  /*!< , */
#define HALFPI      (PI / 2.0)                                  /*!< , */
#define TWOPI       (2.0 * PI)                                  /*!< , */
#define SQRT2       1.4142135623730950488016887242096980785697  /*!< , */
#define SQRT3       1.7320508075688772935274463415058723669428  /*!< , */

/*      Angle conversions */
#define DEG2RAD     (PI / 180.0)        /*!< degrees to radians */
#define RAD2DEG     (180.0 / PI)        /*!< radians to degrees */

/*!     A very small number, used to avoid divide by 0 errors */
#define SFA         1E-10

/*------------------------------------------------------------------------------
 *
 * Other useful macros and inline functions
 *
 *----------------------------------------------------------------------------*/
/*      Convert angle from one unit to another */
#ifdef PREDEF_STANDARD_C_1999
/*          Compiler supports inline functions */
/*! Returns \a angle_deg converted from degrees to radians */
static inline double degToRad(double angle_deg) { return angle_deg * DEG2RAD; }
/*! Returns \a angle_rad converted from radians to degrees */
static inline double radToDeg(double angle_rad) { return angle_rad * RAD2DEG; }

#else
/*          C89/C90 compiler only - no inline functions. Need macros instead */
#define degToRad(angle_deg__)       ((angle_deg__) * DEG2RAD)
#define radToDeg(angle_rad__)       ((angle_rad__) * RAD2DEG)

#endif

/*! Because C passes arrays to functions by passing only a pointer to the
 *  zero'th element of the array, we often need to pass the array size in a
 *  separate argument to the function. Make that easier. */
#define ARRAY_SIZE(x__)             (int)(sizeof(x__)/sizeof(x__[0]))

/*! There are times when we want to declare arrays or buffers whose dimension is
 *  a power of 2, so that the index value can be bit-masked to enforce wrap-
 *  around. Code based on this assumption will fail if the array size is changed
 *  by someone who doesn't notice this requirement.
 *  Here is a check you can use to help prevent that happening */
#define ISPOWER2(x__)               (!((x__)&((x__)-1)))

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/*!     Useful macros for supporting a "Design by Contract" approach to 
        programming embedded systems. These are heavily influenced by the
        article "Design by Contract (DbC) for Embedded Software" by Miro Samek,
        to be found on the Barr Group website, at
        https://barrgroup.com/Embedded-Systems/How-To/Design-by-Contract-for-Embedded-Software

        To use, 
        1. define an implementation of onAssert__() function as specified below.
           For example, it could write the first n chars of the filename and the
           line number to a bit of memory that won't be erased on a reboot, and
           then force a reboot. 
        2. spread REQUIRE(), ENSURE(), and INVARIANT() macros throughout your
           code. (In brief, REQUIRE() macros are generally put at the beginning
           of functions to check the validity of passed parameters. ENSURE()
           macros are generally placed at the end of functions to check that the
           functions have not calculated invalid results.)
        3. put the DEFINE_THIS_FILE statement near the top of each C file that
           contains any REQUIRE(), ENSURE(), or INVARIANT() macros.
        4. choose the ASSERTION_LEVEL that you want.

        (If you define USE_STANDARD_ASSERT as described below, you don't really
        need to do steps 1 and 3. They are really for embedded systems.)
*/
//#define ASSERTION_LEVEL     0       /* DbC assertions are all disabled */
//#define ASSERTION_LEVEL     1       /* Pre-conditions checked (only) */
#define ASSERTION_LEVEL     2       /* Pre & post-conditions checked */
//#define ASSERTION_LEVEL     3       /* Everything checked */

/* Uncomment the following line if you want to use the standard assert() macro,
 * and NOT use your own onAssert__() function to implement these DbC macros.
 * Comment it out if you want to create an onAssert__() function, or if assert.h
 * is not available.  */
#define USE_STANDARD_ASSERT

#if (ASSERTION_LEVEL == 0)
/*!    Leave DEFINE_THIS_FILE undefined. But we get a -Wextra-semi warning
 *     message if it is defined as nothing at all, so use this dummy
 *     declaration instead. */
#  define DEFINE_THIS_FILE    typedef int dummy_t
#  define ASSERT(ignore_)     ((void)0)     /*!< all assertions disabled */

#else

#  if defined(USE_STANDARD_ASSERT)
#    include <assert.h>
     /*!    Leave DEFINE_THIS_FILE undefined. But we get a -Wextra-semi warning
      *     message if it is defined as nothing at all, so use this dummy
      *     declaration instead. */
#    define DEFINE_THIS_FILE    typedef int dummy_t
#    define ASSERT(test_)       assert(test_)   /*!< Uses standard assert */
#  else
#    ifdef __cplusplus
     extern "C"
     {
#    endif
         /*! callback invoked in case of assertion failure */
         void onAssert__(const char *file, unsigned line);
#    ifdef __cplusplus
     }
#    endif

/*!       Ensure exe image contains one copy of file name, not one per ASSERT */
#    define DEFINE_THIS_FILE    static char const THIS_FILE__[] = __FILE__
#    define ASSERT(test_)       ((test_) ? (void)0  \
                                         : onAssert__(THIS_FILE__, __LINE__))
#  endif
#endif

#if (ASSERTION_LEVEL >= 1)
#  define REQUIRE(test_)      ASSERT(test_) /*!< Check preconditions */
#else
#  define REQUIRE(test_)      ((void)0)     /*!< Pre-condition checks disabled*/ 
#endif

#if (ASSERTION_LEVEL >= 2)
#  define ENSURE(test_)       ASSERT(test_) /*!< Check post-conditions */
#else
#  define ENSURE(test_)       ((void)0)     /*!<Post-condition checks disabled*/
#endif

#if (ASSERTION_LEVEL == 3)
#  define INVARIANT(test_)    ASSERT(test_) /*!< Check invariants */
#else
#  define INVARIANT(test_)    ((void)0)     /*!< Invariant checks disabled */
#endif

#if !defined(ASSERTION_LEVEL) || (ASSERTION_LEVEL < 0) || (ASSERTION_LEVEL > 3)
#  error "ASSERTION_LEVEL must be set to one of 0, 1, 2, or 3"
#endif

/* Finally, one more piece of assertion checking that is particularly relevant
 * to C code. Unlike C++, C does not allow passing parameters to functions by
 * reference. It allows passing by address/by pointer, but this is not exactly
 * the same thing. It is quite possible to pass NULL to any parameter in C that
 * is expecting an address or pointer, and the compiler will not object. After
 * all, sometimes functions check for NULL pointers and take different actions
 * if one is passed. But all too often, functions cannot handle getting a NULL
 * pointer, and will crash if one is passed. An assertion check helps here.
 *      Uncomment the following line if you want to perform NULL checking.
 * Comment it out if you want to suppress NULL checking. */
#define ENABLE_NULL_CHECKING

#ifdef ENABLE_NULL_CHECKING
#  include <stddef.h>                       /* Make sure NULL is defined */
#  define REQUIRE_NOT_NULL(pointer_)  ASSERT((pointer_) != NULL)    /*!< , */
#else
#  define REQUIRE_NOT_NULL(pointer_)  ((void)0) /*!< Pointer checking disabled*/
#endif

#endif /*GENERAL_H*/
