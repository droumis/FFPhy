/* $Id: CSStorage.cc,v 1.1 2008/02/01 19:46:49 chengs Exp $

   Sen Cheng, Tue Jun  6 14:06:42 PDT 2006
   definition of class CSStorage
*/

#include "CSStorage.h"

using namespace AFilter;

void CSStorage::init(int start, int end, CSFunction *f) {
    super= f;
    startindex= start; endindex= end;
    T= endindex+1;
    nFct= f->getNFct();
    nPerFct= f->getNPerFct();
    nUpdate= f->getNUpdate();
    nTotal= f->getNTotal();
    const vector<CSFct *> *tmp= super->getFunctionObjs();
    param= new double*[nFct];
    for(int i=0; i<nFct; i++) {
        param[i]= tmp->at(i)->param;
    }
}

