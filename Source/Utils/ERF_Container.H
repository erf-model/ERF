#ifndef ERF_CONTAINER_H
#define ERF_CONTAINER_H

template<typename V, typename T>
bool containerHasElement (const V& iterable, const T& query) {
    return std::find(iterable.begin(), iterable.end(), query) != iterable.end();
}

#endif
