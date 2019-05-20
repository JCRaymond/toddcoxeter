#include <initializer_list>
#include <iostream>
#include <vector>
#include <string>

struct relation {
   int len;
   int* terms;
   relation(): len(0), terms(nullptr) {}
   relation(int* terms, int size) {
      initialize(terms, size);
   }
   relation(std::initializer_list<int> il_terms) {
      int* terms = new int[il_terms.size()];
      int i = 0;
      for (auto it = std::begin(il_terms); it != std::end(il_terms); ++i, ++it) {
         terms[i] = *it;
      }
      initialize(terms, il_terms.size());
      delete[] terms;
   }
   private:

   void initialize(int* terms, int size) {
      int mult = 1;
      int last = terms[size-1];
      if (last < 0) {
         mult = -last;
         size--;
      }
      this->len = size*mult;
      this->terms = new int[this->len];
      int i = 0;
      for (int j = 0; j < mult; j++) {
         for (int k = 0; k < size; k++, i++) {
            this->terms[i] = terms[k];
         }
      }
   }
};

struct group {
   int ngens;
   int nrels;
   relation* rels;
   group(int ngens): ngens(ngens), nrels(0), rels(nullptr) {}
   group(int ngens, std::initializer_list<relation> rels): ngens(ngens), nrels(rels.size()) {
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

struct path_elem {
   int from;
   int gen;
   path_elem(): from(-1), gen(-1) {}
   path_elem(int from, int gen): from(from), gen(gen) {}
};

struct word {
   int* chars;
   int len;
   word(): chars(nullptr), len(0) {}
   word(int* chars, int len): chars(chars), len(len) {}
   word append(word other) {
      int newlen = this->len + other.len;
      int* newchars = new int[newlen];
      int p = 0;
      for (int i = 0; i < this->len; i++) {
         newchars[p++] = this->chars[i];
      }
      for (int i = 0; i < other.len; i++) {
         newchars[p++] = other.chars[i];
      }
   }
};

#define DEFAULT_ALLOC_SIZE 32

struct cosets {
   std::vector<path_elem> path;
   std::vector<int>* gens;
   int ngens;
   int num_cosets;
   int next_ent_cs;
   int next_ent_gn;
   int alloc_size;
   cosets(int ngens): ngens(ngens), num_cosets(1), next_ent_cs(0), next_ent_gn(0), alloc_size(DEFAULT_ALLOC_SIZE) {
      this->path = std::vector<path_elem>(this->alloc_size);
      this->path[0] = path_elem(0, -1);
      if (ngens == 0)
         this->gens = nullptr;
      else
         this->gens = new std::vector<int>[ngens];
      for (int i=0; i<ngens; i++) {
         this->gens[i] = std::vector<int>(this->alloc_size,-1);
      }
   }
   ~cosets() {
      if (this->gens != nullptr)
         delete[] this->gens;
   }
   int apply(int coset, int gen) {
      return this->gens[gen][coset];
   }
   int apply(int coset, word w) {
      for (int i = 0; i < w.len; i++) {
         coset = this->apply(coset, w.chars[i]);
      }
      return coset;
   }
   int apply(int coset, word w, int* gen_map) {
      for (int i = 0; i < w.len; i++) {
         coset = this->apply(coset, gen_map[w.chars[i]]);
      }
      return coset;
   }
   void set(int coset, int gen, int target) {
      this->gens[gen][coset] = target;
      this->gens[gen][target] = coset;
   }
   void find_next_entry() {
      while (this->next_ent_cs < this->num_cosets and this->apply(this->next_ent_cs, this->next_ent_gn) > -1) {
         this->next_ent_gn++;
         if (this->next_ent_gn >= this->ngens) {
            this->next_ent_gn %= this->ngens;
            this->next_ent_cs++;
         }
      }
   }
   bool add_coset() {
      this->find_next_entry();
      if (this->next_ent_cs >= this->num_cosets) {
         this->path.resize(this->num_cosets);
         for (int gen = 0; gen < this->ngens; gen++) {
            this->gens[gen].resize(this->num_cosets);
         }
         return false;
      }
      int next_coset = this->num_cosets;
      this->num_cosets++;
      if (this->alloc_size < this->num_cosets) {
         this->alloc_size <<= 1;
         for (int gen = 0; gen < this->ngens; gen++) {
            this->gens[gen].resize(this->alloc_size, -1);
         }
         this->path.resize(this->alloc_size);
      }
      this->path[next_coset] = path_elem(this->next_ent_cs, this->next_ent_gn);
      this->set(this->next_ent_cs, this->next_ent_gn, next_coset);
      return true;
   }
   word* path_to_words(int* gen_map = nullptr) {
      word* words = new word[this->num_cosets];
      for (int i = 0; i < this->num_cosets; i++) {
         words[i] = word();
      }
      int size, curr;
      for (int i = this->num_cosets - 1; i > 0; i--) {
         if (words[i].len != 0)
            continue;
         size = 0;
         curr = i;
         while (curr != 0) {
            size++;
            curr = this->path[curr].from;
         }
         int* chars = new int[size];
         size--;
         curr = i;
         if (gen_map == nullptr) {
            while (curr != 0) {
               words[curr] = word(chars, size+1);
               chars[size--] = this->path[curr].gen;
               curr = this->path[curr].from;
            }
         }
         else {
            while (curr != 0) {
               words[curr] = word(chars, size+1);
               chars[size--] = gen_map[this->path[curr].gen];
               curr = this->path[curr].from;
            } 
         }
      }
      return words;
   }
   void print(int* gen_map = nullptr) {
      std::cout<<"    ";
      for (int j = 0; j < this->ngens; j++) {
         if (gen_map == nullptr)
            std::cout<<j<<' ';
         else
            std::cout<<gen_map[j]<<' ';
      }
      std::cout<<std::endl;
      for (int i = 0; i < this->num_cosets; i++) {
         std::cout<<i<<" | ";
         for (int j = 0; j < this->ngens; j++) {
            std::cout<<this->apply(i,j)<<' ';
         }
         std::cout<<std::endl;
      }
   }
};

struct reltable_row {
   int left_coset;
   int right_target;
   int* left_gen;
   int* right_gen;
   reltable_row* next;
   reltable_row(): left_coset(0), right_target(0), left_gen(nullptr), right_gen(nullptr), next(nullptr) {}
   reltable_row(relation* rel, int i): next(nullptr) {
      this->left_coset = i;
      this->right_target = i;
      this->left_gen = &(rel->terms[0]);
      this->right_gen = &(rel->terms[rel->len-1]);
   }
   bool apply_learn(cosets* c) {
      while (this->left_gen < this->right_gen) {
         int left_target = c->apply(this->left_coset, *this->left_gen);
         if (left_target < 0)
            break;
         this->left_gen++;
         this->left_coset = left_target;
      }
      while (this->left_gen < this->right_gen) {
         int right_coset = c->apply(this->right_target, *this->right_gen);
         if (right_coset < 0)
            break;
         this->right_gen--;
         this->right_target = right_coset;
      }
      if (this->right_gen == this->left_gen) {
         c->set(this->left_coset, *this->left_gen, this->right_target);
         return true;
      }
      return false;
   }
};

struct reltable {
   group* g;
   reltable_row* first;
   reltable_row* last;
   reltable_row* unused;
   int nrows;
   reltable(group* g): g(g), first(nullptr), last(nullptr), unused(nullptr), nrows(0) {}
   ~reltable() {
      reltable_row* curr, *next;
      curr = this->first;
      while (curr != nullptr) {
         next = curr->next;
         delete curr;
         curr = next;
      }
      curr = this->unused;
      while (curr != nullptr) {
         next = curr->next;
         delete curr;
         curr = next;
      }
   }
   void add_rows() {
      if (this->g->nrels == 0)
         return;
      int rownum = 0;
      reltable_row* row;
      if (this->unused != nullptr) {
         row = this->unused;
         this->unused = this->unused->next;
         *row = reltable_row(&(this->g->rels[0]), this->nrows);
      }
      else {
         row = new reltable_row(&(this->g->rels[0]), this->nrows);
      }
      rownum++;
      if (this->first == nullptr)
         this->first = row;
      else
         this->last->next = row;
      while (this->unused != nullptr and rownum < this->g->nrels) {
         row->next = this->unused;
         this->unused = this->unused->next;
         *(row->next) = reltable_row(&(this->g->rels[rownum]), this->nrows);
         row = row->next;
         rownum++;
      }
      for (; rownum < this->g->nrels; rownum++) {
         row->next = new reltable_row(&(this->g->rels[rownum]), this->nrows);
         row = row->next;
      }
      this->last = row;
      this->nrows++;
   }
   bool apply_learn(cosets* c) {
      bool learned = false;
      reltable_row* next;
      while (this->first != nullptr and this->first->apply_learn(c)) {
         learned = true;
         next = this->first->next;
         this->first->next = this->unused;
         this->unused = this->first;
         this->first = next;
      }
      if (this->first == nullptr) {
         this->last = nullptr;
         return learned;
      }
      reltable_row* curr = this->first;
      while (curr->next != nullptr) {
         next = curr->next;
         if (next->apply_learn(c)) {
            learned = true;
            curr->next = next->next;
            next->next = this->unused;
            this->unused = next;
         }
         else
            curr = curr->next;
      }
      this->last = curr;
      return learned;
   }
};

struct sub_coxeter;

struct coxeter : group {
   coxeter(int ngens): group(ngens) {}
   coxeter(int ngens, relation* rels, int size): group(ngens) {
      initialize(rels, size);
   }
   coxeter(int ngens, std::initializer_list<relation> il_rels): group(ngens) {
      relation* rels = new relation[il_rels.size()];
      int i = 0;
      for (auto it = std::begin(il_rels); it != std::end(il_rels); ++i, ++it) {
         rels[i] = *it;
      }
      initialize(rels, il_rels.size());
      delete[] rels;
   }
   cosets* enumerate_cosets(std::initializer_list<int> gens) {
      if (this->ngens == 0)
         return new cosets(this->ngens);
      cosets* c = new cosets(this->ngens);
      reltable rts(this);
      for (auto it = std::begin(gens); it != std::end(gens); ++it) {
         c->set(0, *it, 0);
      }

      do {
         rts.add_rows();
         while (rts.apply_learn(c));
      } while(c->add_coset());
      
      return c;
   }
   cosets* enumerate_cosets(int* gens, int n) {
      if (this->ngens == 0)
         return new cosets(this->ngens);
      cosets* c = new cosets(this->ngens);
      reltable rts(this);

      for (int i = 0; i < n; i++) {
         c->set(0, gens[i], 0);
      }

      do {
         rts.add_rows();
         while (rts.apply_learn(c));
      } while(c->add_coset());
      
      return c;
   }
   cosets** enumerate_all_cosets() {
      int num_subgroups = 1 << this->ngens;
      int temparr[this->ngens];
      int n, mask; 
      cosets** ec = new cosets*[num_subgroups];
      for (int i = 0; i < num_subgroups; i++) {
         n = 0;
         mask = 0;
         for (int j = 0; j < this->ngens; j++) {
            mask = 1 << j;
            if (i & mask) {
               temparr[n++] = j;
            }
         }
         ec[i] = this->enumerate_cosets(temparr, n);
      }
      return ec;
   }
   sub_coxeter* get_subgroup(int* gens, int m);
   private:

   void initialize(relation* rels, int size) {
      this->nrels = ( this->ngens * (this->ngens - 1) ) / 2;
      if (this->nrels == 0) {
         this->rels = nullptr;
         return;
      }
      this->rels = new relation[this->nrels];
      bool* hasrel[this->ngens];
      for (int i = 0; i < this->ngens; i++) {
         hasrel[i] = new bool[i+1];
         for (int j = 0; j <= i; j++) {
            hasrel[i][j] = false;
         }
      }
      int p = 0;
      relation r;
      for (int i = 0; i < size; i++) {
         r = rels[i];
         int f = r.terms[0];
         int t = r.terms[1];
         if (t == f)
            continue;
         this->rels[p++] = r;
         if (f < t)
            std::swap(f,t);
         hasrel[f][t] = true;
      }
      for (int i = 0; i < this->ngens; i++) {
         for (int j = 0; j <= i; j++) {
            if (i != j and !hasrel[i][j])
               this->rels[p++] = {j, i, -2};
         }
      }
      for (int i = 0; i < this->ngens; i++) {
         delete[] hasrel[i];
      }
   }
};

struct sub_coxeter : coxeter {
   int* gen_map;
   sub_coxeter(int ngens): coxeter(ngens) {}
   sub_coxeter(int ngens, relation* rels, int size, int* gen_map): coxeter(ngens, rels, size), gen_map(gen_map) {}
};

sub_coxeter* coxeter::get_subgroup(int* gens, int m) {
   int nrels = (m  * (m - 1) ) / 2;
   relation* rels = new relation[nrels];
   int gen_map[this->ngens];
   for (int i = 0; i < this->ngens; i++) {
      gen_map[i] = -1;
   }
   for (int i = 0; i < m; i++) {
      gen_map[gens[i]] = i;
   }
   int p = 0;
   relation r;
   int f,t;
   for (int i = 0; i < this->nrels; i++) {
      r = this->rels[i];
      f = gen_map[r.terms[0]];
      t = gen_map[r.terms[1]];
      if (f != -1 and t != -1) {
         rels[p++] = {f, t, -r.len/2};
      }
   }
   return new sub_coxeter(m, rels, nrels, gens);
}

void print_words(word* words, int n, int* gen_map = nullptr) {
   word w;
   for (int i = 0; i < n; i++) {
      w = words[i];
      std::cout<<i<<": ";
      for (int j = 0; j < w.len; j++) {
         if (gen_map != nullptr)
            std::cout<<gen_map[w.chars[j]];
         else
            std::cout<<w.chars[j];
      }
      std::cout<<std::endl;
   }
}

struct coset_points {
   int* mapping;
   int id;
   coset_points(): mapping(nullptr), id(-1) {}
   coset_points(int* mapping, int id): mapping(mapping), id(id) {}
};

std::vector<coset_points>* thing(coxeter* c){
   auto points = c->enumerate_cosets({});
   int num_subgroups = 1 << c->ngens;
   
   auto things = new std::vector<coset_points>[c->ngens + 1];
   int gens[c->ngens];
   int n, mask;
   for (int i = 0; i < num_subgroups; i++) {
      n = 0;
      mask = 1;
      for (int j = 0; j < c->ngens; j++) {
         if (i & mask) {
            gens[n++] = j;
         }
         mask <<= 1;
      }

      auto subgroup = c->get_subgroup(&gens[0], n);
      
      auto subgroup_cosets = c->enumerate_cosets(gens, n);
      auto subgroup_coset_words = subgroup_cosets->path_to_words();

      auto subgroup_points = subgroup->enumerate_cosets({});
      auto subgroup_words = subgroup_points->path_to_words(subgroup->gen_map);

      int* group_mapping = new int[points->num_cosets];
      int temppoint;
      for (int point = 0; point < subgroup_points->num_cosets; point++) {
         temppoint = points->apply(0, subgroup_words[point]);
         for (int coset = 0; coset < subgroup_cosets->num_cosets; coset++) {
            group_mapping[points->apply(temppoint, subgroup_coset_words[coset])] = coset;
         }
      }
      things[n].push_back(coset_points(group_mapping, i));
   }
   return things;
}

struct primitive {
   int dim;
   int* p;
   primitive(int dim): dim(dim), p(new int[dim+1]) {}
   primitive(int dim, int* arr): primitive(dim) {
      for (int i = 0; i <= dim; i++) {
         p[i] = arr[i];
      }
   }
   primitive(std::initializer_list<int> vals): primitive(vals.size()-1) {
      int i = 0;
      for (auto it = std::begin(vals); it != std::end(vals); ++it) {
         p[i++] = *it;
      }
   }
   primitive* copy() {
      return new primitive(dim, p);
   }
   primitive* reorient_() {
      if (dim > 0)
         std::swap(p[0], p[1]);
      return this;
   }
   primitive* reorient() {
      return this->copy()->reorient();
   }
   primitive* shift_(cosets* c, word w) {
      for (int i = 0; i <= dim; i++) 
         p[i] = c->apply(p[i], w);
      return this;
   }
   primitive* shift(cosets* c, word w) {
      return this->copy()->shift_(c, w);
   }
   primitive* grow(int elem) {
      primitive* n = new primitive(dim + 1, p);
      n->set(dim+1,elem);
      return n;
   }
   primitive* shift_grow(cosets* c, word w, int elem) {
      return this->grow(0)->shift_(c, w)->set(dim+1, elem);
   }
   ~primitive() {
      delete[] p;
   }
   int get(int pos) {
      return p[pos];
   }
   primitive* set(int pos, int val) {
      p[pos] = val;
      return this;
   }
   void print() {
      std::cout<<'(';
      for (int i = 0; i <= dim; i++) {
         std::cout<<p[i];
         if (i != dim)
            std::cout<<", ";
      }
      std::cout<<')'<<std::endl;
   }
};

int choose(int n, int k) {
   if (k > n - k) {
      k = n - k;
   }
   int a_;
   int res = 1;
   int b = 2;
   for (int a = n; a > n - k; a--) {
      a_ = a;
      while (b <= k and a % b == 0) {
         a_ /= b;
         b++;
      }
      res *= a_;
      while (b <= k and res % b == 0) {
         res /= b;
         b++;
      }
   }
   return res;
}

int main(){
   //coxeter<6> g = {{0,1,-3},{1,2,-3},{2,3,-3},{2,4,-3},{4,5,-3}};
   //coxeter<7> g = {{0,1,-3},{1,2,-3},{2,3,-3},{2,4,-3},{4,5,-3},{5,6,-3}};
   //coxeter g(4,{{0,1,-5},{1,2,-3},{2,3,-3}});
   //coxeter g(3,{{0,1,-4},{1,2,-3}});
   //coxeter g(3, {{0,1,-3}});
   coxeter g(2, {{0,1,-4}});
   cosets* t_c = g.enumerate_cosets({1});

   int num_subgroups = 1 << g.ngens;

   int indexes[g.ngens + 1];
   indexes[0] = 0;
   for (int i = 1; i <= g.ngens; i++) {
      indexes[i] = indexes[i-1] + choose(g.ngens, i-1);
   }

   int ids[num_subgroups];
   int dims[num_subgroups];
   int* sub_gens[num_subgroups];
   int tmpgens[g.ngens];
   int n, mask, index, *gens;
   for (int i = 0; i < num_subgroups; i++) {
      n = 0;
      mask = 1;
      for (int j = 0; j < g.ngens; j++) {
         if (i & mask)
            tmpgens[n++] = j;
         mask <<= 1;
      }
      gens = new int[n];
      for (int j = 0; j < n; j++) {
         gens[j] = tmpgens[j];
      }
      index = indexes[n];
      indexes[n]++;
      ids[index] = i;
      dims[index] = n;
      sub_gens[index] = gens;
   }
   
   std::vector<primitive*>* primitives = new std::vector<primitive*>[num_subgroups];
   std::vector<primitive*>* final_primitives = new std::vector<primitive*>[g.ngens + 1];

   primitive *p;
   p = new primitive({0});
   primitives[0].push_back(p);
   final_primitives[0].push_back(p);

   std::vector<primitive*> sg_primitives, *curr_primitives, *dim_primitives;
   int sgens[g.ngens];
   int id, sub_id, m, si1, si2, num_prims;
   bool ro_si1;
   sub_coxeter *sg;
   cosets *c;
   word *c_w, w;
   for (int i = 0; i < num_subgroups; i++) {
      id = ids[i];
      n = dims[i];
      gens = sub_gens[i];
      m = n - 1;
      
      sg = g.get_subgroup(gens, n);

      curr_primitives = &(primitives[id]);
      dim_primitives = &(final_primitives[n]);

      if (n == 2) {
         si1 = 1 << gens[0];
         si2 = id & ~si1;
         if (si1 > si2) {
            std::swap(si1, si2);
         }
         ro_si1 = false;
         if ((si2 - si1) % 2 == 1)
            ro_si1 = true;
         
         c = sg->enumerate_cosets({0});
         c_w = c->path_to_words(sg->gen_map);
         sg_primitives = primitives[si1];
         num_prims = sg_primitives.size();
         for (int j = 1; j < c->num_cosets; j++) {
            w = c_w[j];
            for (int k = 0; k < num_prims; k++) {
               p = sg_primitives[k]->shift_grow(t_c, w, 0);
               if (w.len % 2 == 1 xor ro_si1)
                  p->reorient_();
               curr_primitives->push_back(p);
               dim_primitives->push_back(p);
            }
         }

         c = sg->enumerate_cosets({1});
         c_w = c->path_to_words(sg->gen_map);
         sg_primitives = primitives[si2];
         num_prims = sg_primitives.size();
         for (int j = 1; j < c->num_cosets; j++) {
            w = c_w[j];
            for (int k = 0; k < num_prims; k++) {
               p = sg_primitives[k]->shift_grow(t_c, w, 0);
               if (w.len % 2 == 1 xor !ro_si1)
                  p->reorient_();
               curr_primitives->push_back(p);
               dim_primitives->push_back(p);
            }
         }

      }
      else {
         for (int j = 0; j < n; j++)
            sgens[j] = j;
         for (int j = m; j >= 0; j--) {
            std::swap(sgens[j], sgens[m]);
            sub_id = id & ~(1 << gens[j]);
            
            c = sg->enumerate_cosets(sgens, m);
            c_w = c->path_to_words(sg->gen_map);
            sg_primitives = primitives[sub_id];
            num_prims = sg_primitives.size();
            for (int k = 1; k < c->num_cosets; k++) {
               w = c_w[k];
               for (int l = 0; l < num_prims; l++) {
                  p = sg_primitives[l]->shift_grow(t_c, w, 0);
                  if (w.len % 2 == 1)
                     p->reorient_();
                  curr_primitives->push_back(p);
                  dim_primitives->push_back(p);
               }
            }

         }
      }
      
      /*
      std::cout<<id<<':'<<std::endl;
      for (int j = 0; j < curr_primitives->size(); j++) {
         std::cout<<'\t';
         (*curr_primitives)[j]->print();
      }
      std::cout<<std::endl;
      */
      
      c = g.enumerate_cosets(gens, n);
      c_w = c->path_to_words();
      sg_primitives = *curr_primitives;
      num_prims = curr_primitives->size();
      for (int j = 1; j < c->num_cosets; j++) {
         w = c_w[j];
         for (int k = 0; k < num_prims; k++) {
            p = sg_primitives[k]->shift(t_c, w);
            if (w.len % 2 == 1)
               p->reorient_();
            dim_primitives->push_back(p);
         }
      }
   }

   for (int i = 0; i <= g.ngens; i++) {
      sg_primitives = final_primitives[i];
      std::cout<<i<<':'<<std::endl;
      for (int j = 0; j < sg_primitives.size(); j++) {
         std::cout<<'\t';
         sg_primitives[j]->print();
      }
      std::cout<<std::endl;
   }

   /*
   std::vector<coset_points> temp;
   int* thingy;
   for (int i = 0; i < g.ngens + 1; i++) {
      temp = res[i];
      std::cout<<i<<" dim: "<<std::endl;
      for (int j = 0; j < temp.size(); j++) {
         thingy = temp[j].mapping;
         std::cout<<temp[j].id<<": ";
         for (int k = 0; k < points->num_cosets; k++) {
            std::cout<<thingy[k]<<' ';
         }
         std::cout<<std::endl;
      }
   }
   */
   return 0;
}
