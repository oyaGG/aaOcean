//Copyright(C) 2006, William Chan
//All rights reserved.
//
//    Redistribution and use in source and binary forms, with or without
//  modification, are permitted provided that the following conditions are met:
//  
//    1) Redistributions of source code must retain the above copyright
//      notice, this list of conditions and the following disclaimer.
//    2) Redistributions in binary form must reproduce the above copyright
//      notice, this list of conditions and the following disclaimer in the
//      documentation and/or other materials provided with the distribution.
//    3) Redistributions of source code must be provided at free of charge.
//    4) Redistributions in binary forms must be provided at free of charge.
//    5) Redistributions of source code within another distribution must be
//      provided at free of charge including the distribution which is
//      redistributing the source code. Also, the distribution which is
//      redistributing the source code must have its source code
//      redistributed as well.
//    6) Redistribution of binary forms within another distribution must be
//      provided at free of charge.
//      
//
//    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
//  IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
//  THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
//  PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS
//  BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
//  CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
//  SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
//  INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
//  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
//  ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
//  POSSIBILITY OF SUCH DAMAGE.
//
//-------------------------------------------------------------------------------
//
//Website: http://williamchan.ca
//
//I would appreciate it if you email me at me@williamchan.ca if you
//use my source code in any of your projects. Please do not hesitate to leave me
//comments and suggestions too.
#ifndef ALIGNED_MALLOC
#define ALIGNED_MALLOC
#include "stdlib.h"
 
void* aligned_malloc(size_t bytes, size_t alignment)
{
    // call standard c malloc to allocate memory
    // aloocate memory size = bytes + alignment + sizeof(alignment)
    // bytes               - obvious to allocate at least the user requested amount
    // alignment           - this is the padding in front or behind or both the memory block
    // sizeof(alignment)   - this is to save the real value pointer value (ptr vs aligned_ptr) before the aligned memory block
    //                     - this is necessary to free the entire memory block and not just the aligned memory block
    //                     - sizeof(size_t) was chosen over sizeof(unsigned int) for compatibility of 32 and 64-bit architectures
	
	if(alignment == 0 || alignment == 1)			//1 and 0 means no alignment needed
		return malloc(bytes);						//regular malloc
	else if(((alignment - 1) & alignment) != 0)		//check if power of 2
		return NULL;								//not a power of 2...failure
 
	uintptr_t ptr = (uintptr_t)malloc(bytes + alignment + sizeof(uintptr_t)); //allocate enough memory for everything
 
	if(ptr == NULL) //failed malloc! :(
		return NULL;
 
	uintptr_t aligned_ptr = (ptr + alignment + sizeof(uintptr_t)) & (~(alignment - 1)); //align!
	*(uintptr_t*)(aligned_ptr - sizeof(uintptr_t)) = ptr; //store the entire memory block address before the aligned block
 
	return (void*)aligned_ptr;
}
 
void aligned_free(void* p)
{
	free((void*)(*((uintptr_t*)((uintptr_t)p - sizeof(uintptr_t))))); //get the pointer of the entire block and then FREE!
}

#endif /* ALIGNED_MALLOC_H */