#include "klu_internal.h"

void push_back(list* l, int v){
    if(l->length==0){
        node* n = (node*)malloc(sizeof(node));
        n->value = v;
        n->next = NULL;
        n->prev = NULL;
        // first push
        l->head = n;
        l->tail = n;
    } else {
        node* n = (node*)malloc(sizeof(node));
        n->value = v;
        n->next = NULL;
        n->prev = NULL;
        node* t = l->tail;
        l->tail = n;
        n->prev = t;
        t->next = n;
    }
    l->length++;
    return;
}

void push_front(list* l, int v){
    node* n = (node*)malloc(sizeof(node));
    n->value = v;
    n->next = NULL;
    n->prev = NULL;
    if(l->length==0){
        // first push
        l->head = n;
        l->tail = n;
    } else {
        node* h = l->head;
        l->head = n;
        n->prev = h;
        h->next = n;
    }
    l->length++;
    return;
}

int pop_front(list* l){
    if(l->length==0||!l->head)
        return 0;
    node* h = l->head;
    node* n = h->next;
    int ret = h->value;
    h->next = NULL;
    h->prev = NULL;
    if(n){
        l->head = n;
        n->prev = NULL;
    } else {
        // h has no next -> only node in list
        l->head = NULL;
        l->tail = NULL;
    }
    l->length--;
    return ret;
}

int pop_back(list* l){
    if(l->length==0||!l->tail)
        return 0;
    node* t = l->tail;
    node* p = t->prev;
    int ret = t->value;
    t->next = NULL;
    t->prev = NULL;
    if(p){
        l->tail = p;
        p->next = NULL;
    } else {
        l->tail = NULL;
        l->head = NULL;
    }
    l->length--;
    return ret;
}

void concat(list* left, list* right){
    node* l_t = left->tail;
    node* r_h = right->head;
    if(l_t && r_h){
        l_t->next = r_h;
        r_h->prev = l_t;
        left->tail = right->tail;
        left->length += right->length;
        //unlist(right);
        return;
    } else if(l_t){
        // r_h == NULL
        // do nothing
        return;
    } else if(r_h){
        // l_t == NULL
        // left becomes right
        left->head = right->head;
        left->tail = right->tail;
        left->length = right->length;
        //unlist(right);
        return;
    } else {
        // both are null
        // do nothing
        return;
    }
}

void concat_rev(list* left, list* right){
    node* l_h = left->head;
    node* r_t = right->tail;
    if(l_h && r_t){
        l_h->prev = r_t;
        r_t->next = l_h;
        left->head = right->head;
        left->length += right->length;
        unlist(right);
        return;
    } else if(l_h) {
        return;
    } else if(r_t) {
        left->head = right->head;
        left->tail = right->tail;
        left->length = right->length;
        unlist(right);
        return;
    } else {
        return;
    }
}

/*! \brief remove list, keep nodes
*/
void unlist(list* l){
    l->head = NULL;
    l->tail = NULL;
    l->length = 0;
    free(l);
    return;
}

/*! \brief remove list, remove nodes
*/
void killlist(list* l){
    if(!l||!l->head)
        return;
    node* h = l->head;
    node* tmp;
    while(h){
        tmp = h->next;
        h->next = NULL;
        h->prev = NULL;
        h->value = 0;
        //h = NULL;
        free(h);
        h = tmp;
    }
    return;
}

void printlist(list* l){
    if(l->length==0)
        return;
    node* n = l->head;
    while(n&&n!=l->tail){
        printf("%d->",n->value);
        n = n->next;
    }
    if(n)
        printf("%d",n->value);
    printf("\n");
    return;
}

void printlistrev(list* l){
    if(l->length==0)
        return;
    node* n = l->tail;
    while(n&&n!=l->head){
        printf("%d<-", n->value);
        n = n->prev;
    }
    if(n)
        printf("%d",n->value);
    printf("\n");
    return;
}

void bubblesort(list* l){
    if(!l->head)
        return;
    int piv = pop_front(l);
    quicksort(l,piv);
}

void quicksort(list* l, int piv){
    if(!l)
        return;
    list* left = (list*)malloc(sizeof(list));
    left->head = NULL;
    left->tail = NULL;
    left->length = 0;
    list* right = (list*)malloc(sizeof(list));
    right->head = NULL;
    right->tail = NULL;
    right->length = 0;
    node* iter = l->head;
    node* tmp;
    int iterval = 0;
    while(iter){
        tmp = iter->next;
        iterval = pop_front(l);
        if(iterval < piv){
            push_back(left,iterval);
        } else {
            push_back(right,iterval);
        }
        iter = tmp;
    }
    int front;
    if(left->head){
        front = pop_front(left);
        quicksort(left, front);
    }
    push_back(left, piv);
    if(right->head){
        front = pop_front(right);
        quicksort(right, front);
    }
    concat(left, right);
    concat(l, left);
}

void sort_list(list* l, int method){
    switch(method){
        case BUBBLE: bubblesort(l); break;
        case QUICK: if(l->head) {int init = pop_front(l); quicksort(l, init);} break;
        default: bubblesort(l);
    }
    return;
}


int search(list* l, const int piv){
    node* n = l->head;
    while(n){
        if(n->value == piv){
            return FOUND;
        }
        n = n->next;
    }
    return NOT_FOUND;
}

void push_front_list(list* f, list* b){ /* pushes f in front of b */
    node *back = f->tail;
    node *front = b->head;
    if(back && front){
        back->next = front;
        front->prev = back;
        b->head = f->head;
        f->tail = b->tail;
        f->length += b->length;
        b->length = f->length;
    } else if(back){
        b->head = f->head;
        b->tail = f->tail;
        b->length = f->length;
    } else if(front){
        f->head = b->head;
        f->tail = b->tail;
        f->length = b->length;
    } else {
        ;
    }
}

void push_back_list(list* f, list* b){ /* pushes b back of f */

}