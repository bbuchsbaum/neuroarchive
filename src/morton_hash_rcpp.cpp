#include <Rcpp.h>
#include <openssl/sha.h>
// [[Rcpp::depends(openssl)]]
using namespace Rcpp;

// [[Rcpp::export]]
std::string morton_indices_to_hash_rcpp(IntegerVector voxelIdx) {
    SHA_CTX ctx;
    SHA1_Init(&ctx);
    for (int i = 0; i < voxelIdx.size(); ++i) {
        int val = voxelIdx[i];
        SHA1_Update(&ctx, &val, sizeof(int));
    }
    unsigned char hash[SHA_DIGEST_LENGTH];
    SHA1_Final(hash, &ctx);
    char buf[41];
    for (int i = 0; i < SHA_DIGEST_LENGTH; ++i) {
        std::sprintf(buf + (i*2), "%02x", hash[i]);
    }
    buf[40] = '\0';
    return std::string("sha1:") + std::string(buf);
}
