#include <initializer_list>
#include <iostream>
#include <vector>
#include <iomanip>

using namespace std;

struct relation {
   int len;
   int* terms;
   relation(): len(0), terms(nullptr) {}
   relation(std::initializer_list<int> terms) {
      int size = terms.size();
      int mult = 1;
      int last = *(std::end(terms)-1);
      if (last < 0) {
         mult = -last;
         size--;
      }
      this->len = size*mult;
      this->terms = new int[this->len];
      int i = 0;
      for (int j = 0; j < mult; j++) {
         auto it = std::begin(terms);
         for (int k = 0; k < size; k++, ++it, i++) {
            this->terms[i] = *it;
         }
      }
   }
   void print() {
      for (int i = 0; i < len; i++) {
         cout<<this->terms[i]<<' ';
      }
      cout<<endl;
   }
};

template<int N>
struct group {
   static const int ngens = N;
   int nrels;
   relation* rels;
   group(): nrels(0), rels(nullptr) {}
   group(std::initializer_list<relation> rels): nrels(rels.size()) {
      this->rels = new relation[this->nrels];
      int i = 0;
      for (auto it = std::begin(rels); it != std::end(rels); ++it) {
         this->rels[i++] = *it;
      }
   }
   ~group() {
      delete[] rels;
   }
};

template<int N>
struct coxeter : group<N> {
   coxeter(): group<N>() {}
   coxeter(initializer_list<relation> prels) {
      this->nrels = ( N * (N + 1) ) / 2;
      if (this->nrels == 0) {
         this->rels = nullptr;
         return;
      }
      this->rels = new relation[this->nrels];
      bool* hasrel[N];
      for (int i = 0; i < N; i++) {
         hasrel[i] = new bool[i+1];
         for (int j = 0; j <= i; j++) {
            hasrel[i][j] = false;
         }
      }
      int p = 0;
      for (auto it = std::begin(prels); it != std::end(prels); ++it) {
         relation r = *it;
         int f = r.terms[0];
         int t = r.terms[1];
         if (t == f)
            continue;
         this->rels[p++] = r;
         if (f < t)
            std::swap(f,t);
         hasrel[f][t] = true;
      }
      for (int i = 0; i < N; i++) {
         for (int j = 0; j <= i; j++) {
            if (i == j)
               this->rels[p++] = {i, i};
            else if (!hasrel[i][j])
               this->rels[p++] = {j, i, -2};
         }
      }
      for (int i = 0; i < N; i++) {
         delete[] hasrel[i];
      }
   }
};

struct transform {
   int frw;
   int inv;
   transform(int frw = -1, int inv = -1): frw(frw), inv(inv) {}
};

template<int N>
struct cosets {
   vector<transform>** gens;
   int num_cosets;
   int next_ent_cs;
   int next_ent_gn;
   cosets(): num_cosets(1), next_ent_cs(0), next_ent_gn(0) {
      this->gens = new vector<transform>*[N];
      for (int i=0; i<N; i++) {
         this->gens[i] = new vector<transform>();
      }
   }
   ~cosets() {
      for (int i=0; i<N; i++) {
         delete this->gens[i];
      }
      delete[] this->gens;
   }
   transform* find(int gen, int coset) {
      auto trnsfrms = this->gens[gen];
      if (trnsfrms->size() < this->num_cosets) {
         trnsfrms->resize(this->num_cosets);
      }
      if (this->num_cosets <= coset) {
         return new transform();
      }
      return &(*trnsfrms)[coset];
   }
   transform get(int gen, int coset) {
      return *(this->find(gen, coset));
   }
   int apply(int coset, int gen) {
      return this->find(gen, coset)->frw;
   }
   int inv_apply(int gen, int target) {
      return this->find(gen, target)->inv;
   }
   void find_next_entry() {
      while (this->apply(this->next_ent_cs, this->next_ent_gn) > -1) {
         this->next_ent_gn++;
         if (this->next_ent_gn >= N) {
            this->next_ent_gn %= N;
            this->next_ent_cs++;
         }
      }
   }
   void set(int coset, int gen, int target) {
      this->find(gen, coset)->frw = target;
      this->find(gen, target)->inv = coset;
   }
   bool add_coset() {
      this->find_next_entry();
      if (this->next_ent_cs >= this->num_cosets)
         return false;
      int next_coset = this->num_cosets;
      this->num_cosets++;
      this->set(this->next_ent_cs, this->next_ent_gn, next_coset);
      return true;
   }
   void print() {
      cout<<"    ";
      for (int j = 0; j < N; j++) {
         cout<<j<<' ';
      }
      cout<<endl;
      for (int i = 0; i < this->num_cosets; i++) {
         cout<<i<<" | ";
         for (int j = 0; j < N; j++) {
            cout<<this->apply(i,j)<<' ';
         }
         cout<<endl;
      }
   }
};

template<int N>
int** cosets_to_arr(cosets<N>* c) {
   int** arr = new int*[c->num_cosets];
   for (int i = 0; i < c->num_cosets; i++) {
      int* row = new int[N];
      arr[i] = row;
      for (int j = 0; j < N; j++) {
         row[j] = c->apply(i, j);
      }
   }
   return arr;
}

struct reltable_row {
   relation* rel;
   int* row;
   int rlen;
   int left;
   int right;
   bool complete;
   reltable_row* next;
   reltable_row(): rel(nullptr), row(nullptr), rlen(0), left(0), right(0), complete(true), next(nullptr) {}
   reltable_row(relation* rel, int i): rel(rel), rlen(rel->len+1), left(0), complete(false), next(nullptr) {
      this->row = new int[this->rlen];
      this->row[0] = i;
      this->row[rel->len] = i;
      for (int i = 1; i < rel->len; i++) {
         this->row[i] = -1;
      }
      this->right = rel->len;
   }
   template<int N>
   bool apply_learn(cosets<N>* c) {
      if (this->complete)
         return false;
      while (this->right - this->left > 1) {
         int left_tgt = c->apply(this->row[this->left], this->rel->terms[this->left]);
         if (left_tgt < 0) {
            break;
         }
         this->left++;
         this->row[this->left] = left_tgt;
      }
      while (this->right - this->left > 1) {
         int right_coset = c->inv_apply(this->rel->terms[this->right - 1], this->row[this->right]);
         if (right_coset < 0) {
            break;
         }
         this->right--;
         this->row[this->right] = right_coset;
      }
      if (this->right - this->left == 1) {
         c->set(this->row[this->left], this->rel->terms[this->left], this->row[this->right]);
         this->complete = true;
         return true;
      }
      return false;
   }
   ~reltable_row() {
      if (this->row != nullptr)
         delete[] this->row;
   }
   void print() {
      for (int i = 0; i < this->rlen; i++) {
         cout<<this->row[i]<<' ';
      }
      cout<<endl;
   }
};

struct reltable {
   relation* rel;
   reltable_row* first;
   reltable_row* last;
   reltable_row* comp;
   reltable_row* comp_l;
   int nrows;
   bool complete;
   reltable(relation* rel, int rows = 1): rel(rel), nrows(rows), complete(false) {
      this->first = new reltable_row(this->rel, 0);
      this->last = this->first;
      for (int i = 1; i < rows; i++) {
         this->add_row();
      }
      this->comp = new reltable_row();
      this->comp_l = this->comp;
   }
   void add_row() {
      reltable_row* nrow = new reltable_row(this->rel, this->nrows++);
      if (this->first == nullptr) {
         this->first = nrow;
         this->last = nrow;
         this->complete = false;
      }
      else {
         this->last->next = nrow;
         this->last = nrow;
      }
   }
   template<int N>
   bool apply_learn(cosets<N>* c) {
      if (complete)
         return false;
      bool learned = false;
      while (this->first != nullptr and this->first->apply_learn(c)) {
         this->comp_l->next = this->first;
         this->comp_l = this->first;
         learned = true;
      }
      if (this->first == nullptr) {
         this->last = nullptr;
         this->complete = true;
         return true;
      }
      reltable_row* curr = this->first;
      reltable_row* next;
      while (curr->next != nullptr) {
         next = curr->next;
         if (next->apply_learn(c)) {
            learned = true;
            this->comp_l->next = next;
            this->comp_l = next;
            curr->next = next->next;
         }
         else
            curr = curr->next;
      }
      this->last = curr;
      return learned;
   }
   void print() {
      cout<<' ';
      this->rel->print();
      reltable_row* curr = this->comp->next;
      while (curr != nullptr) {
         curr->print();
         curr = curr->next;
      }
      curr = this->first;
      while (curr != nullptr) {
         curr->print();
         curr = curr->next;
      }
   }
};

template<int N>
cosets<N>* enumerate_cosets(group<N>* g, initializer_list<int> gens) {
   cosets<N>* c = new cosets<N>();
   reltable* rts[g->nrels];
   for (int i = 0; i < g->nrels; i++) {
      rts[i] = new reltable(&(g->rels[i]));
   }

   for (auto it = std::begin(gens); it != std::end(gens); ++it) {
      c->set(0, *it, 0);
   }

   bool learned;
   bool completed;
   do {
      learned = true;
      completed = true;
      while (learned) {
         learned = false;

         for (int i = 0; i < g->nrels; i++) {
            reltable* rt = rts[i];
            if (rt->apply_learn(c))
               learned = true;
            completed = completed and rt->complete;
         }
      }

      if (c->add_coset()) {
         for (int i = 0; i < g->nrels; i++) {
            rts[i]->add_row();
         }
      }
      else
         break;
   } while(!completed);

   return c;
}

int main(){
   coxeter<4> cube = {{0,1,-5},{1,2,-3},{2,3,-3}};
   auto c = enumerate_cosets(&cube,{3});
   c->print(); 
}
