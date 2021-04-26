    
    BITS 64
    
	global _mulmod_64
; callable from C funtion 
; extern "C" uint64_t  _mulmod_64(uint64_t a, uint64_t b, uint64_t m);
; returns (a*b) mod m

;  On Windows, the first four parameters go into rcx, rdx, r8, and r9.
_mulmod_64:
    mov     rax, rcx	; rcx = 1st argument = a
    mov     r11, rdx	; move 2nd arg, rdx=b to r11
	mul		r11			; r11=2nd argument=b,  rdx:rax now holds a*b
    div     r8			; divide rdx:rax=a*b by r8=m, remainder is in rdx, quotient is in rax
    mov     rax, rdx	; move remainder to return register
    ret
