#if (defined(__linux) && !defined(__ia64))
#include <fpu_control.h>
void __attribute__ ((constructor))
     trapfpe () {
  fpu_control_t cw = _FPU_DEFAULT & ~(_FPU_MASK_IM | _FPU_MASK_ZM |
				      _FPU_MASK_OM);
  _FPU_SETCW(cw);
}
void MAIN__()  {}
#endif

void izrtoc_() {}
void igmess_() {}
void igloc2_() {}
void igpxmp_() {}
void izitoc_() {}
/*
void ffinit_() {}
void ffkey_() {}
void ffgo_() {}
*/
void kuproi_() {}
void kuproc_() {}
void kupror_() {}
void kualfa_() {}
void umlog_() {}
void czgeta_() {}
void czputa_() {}
/*======================= hbook dummies ================================
 *
 * From gplmat */
void hdelet_() {}
void hphist_() {}
void hbookb_() {}
void hfill_() {}
void hidopt_() {}
/*
 * From gbhsta */
void hbook1_() {}
void hbookn_() {}
void hcopy_() {}
/*
 * From AliRun */
void hlimit_() {}
/*
 * From HIGZ */

void iacwk_(int* dum) {}
void iclrwk_(int* dum,int* dum2) {}
void iclwk_(int* dum,int* dum2) {}
void idawk_(int* dum) {}
void ifa_(int* n,float* x, float* y) {}
void igbox_(float* x1,float* x2,float* y1,float* y2) {}
void ightor_(float* h,float* l,float* s,
                                    float* r,float* g,float* b) {}
void igpave_(float* x1,float* x2,float* yy1,
                                    float* yy2,float* dum4,int* isbox,
                                    int* isfram,const char* dum5, const int dum) {}
void igpid_(int* dum,const char* name,int* pid,
                                   const char* dum6 , const int l1, const int dum8) {}
void igq_(const char *name,float* rval,  const int l1) {}
void igrng_(float* xsize,float* ysize) {}
void igsa_(int* dum) {}
void igset_(const char *name,float* rval, const int l1) {}
void igterm_() {}
void iopwk_(int* iwkid,int* iconid,int* iwtypi) {}
void ipl_(int* n,float* x,float* y) {}
void ipm_(int* n,float* x,float* y) {}
void irqlc_(int* dum, int* dum2, int* dum3, int* dum4, float* dum5, float* dum6) {}
void iscr_(int* dum,int* ici,float* r,float* g,float* b) {}
void isfaci_(int* col) {}
void isfais_(int* is) {}
void isln_(int* ln) {}
void ismk_(int* mk) {}
void islwsc_(float* wl) {}
void isplci_(int* col) {}
void ispmci_(int* col) {}
void istxci_(int* col) {}
void isvp_(int* dum,float* dum2,float* dum3,float* dum4,float* dum5) {}
void iswn_(int* dum,float* x1,float* x2,float* y1,float* y2) {}
void itx_(float* x,float* y,const char* ptext,  const int l1p) {}
void hplint_(int* dum) {}
void hplend_() {}
void hplfra_(float* x1,float* x2,float* y1, 
                                    float* y2,const char* dum6, const int dum) {}
void igmeta_(int* dum, int* dum2) {}
void iselnt_(int* dum) {}
int igiwty_(int* dum) {return 1;}
void igqwk_(int* dum, const char *name, float* rval, const int l1) {}

