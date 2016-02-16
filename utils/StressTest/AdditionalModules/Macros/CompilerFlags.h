#ifndef COMPILER_FLAGS_H

#define COMPILER_FLAGS_H


// определение используемого компилятора
#if defined(_MSC_VER)

	// признак того, что используется компилятор Microsoft Visual C++
#	define MS_VC_COMPILER

#elif (defined(__GNUC__) || defined(__GNUG__)) && !(defined(__clang__) || defined(__INTEL_COMPILER))

	// признак того, что используется компилятор GNU GCC C/C++
#	define GNU_GCC_COMPILER

#elif defined(__INTEL_COMPILER) || defined(__ICC)

	// признак того, что используется компилятор Intel C/C++
#	define INTEL_C_COMPILER

#elif defined(__clang__)

	// признак того, что используется компилятор Clang
#	define CLANG_COMPILER

#endif


#endif // COMPILER_FLAGS_H