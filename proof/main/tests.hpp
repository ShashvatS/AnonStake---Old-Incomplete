#ifndef TESTS_H
#define TESTS_H

#include "fft.hpp"
#include "link.hpp"

int testfft1(int argc, char *argv[]) {
    init();

    int numProofs = 12;
    if (argc > 1) numProofs = atoi(argv[1]);
    if (numProofs <= 0 || numProofs >= 400) numProofs = 10;

    int numTrials = 1;
    if (argc > 2) numTrials = atoi(argv[2]);
    if (numTrials <= 0 || numTrials >= 100) numTrials = 1;

    // int pause; cin >> pause;

    for (int i = 0; i < numTrials; ++i) {
        cout << "TRIAL # " << i << "\n\n";
        benchmarkProofCreation<FieldR, 120000>(numProofs);
        cout << "\n\n";
    }

    return 0;
}

int testlink1(int argc, char* argv[]) {
    std::cout << "Hello world!\n";
    hello_world();

    std::cout << "Results: " << initParams(5, "testfile.txt") << '\n';
    std::cout << "Results: " << loadParams(false, "testfile.txt") << '\n';

    uint64_t arr[5] = {0, 1, 2, 3, 4};

    std::pair<char*, int> proof = createProofPartial(5, arr);
    std::cout << "Results: " << proof.second << '\n';

    for (int i = 0; i < 32; ++i) {
        std::cout << proof.first[i];
    }
    std::cout << '\n';
    
    //initFieldR();
    // checkFieldR();
    //run();

    FieldR a;
    mp_limb_t *p = (mp_limb_t*)&a;
    *p++ = 1;
    *p++ = 1;
    *p++ = 1;
    *p++ = 1;

    std::cout << sizeof(a) << "\n\n";
    a.print();
    std::cout << a << '\n';

    return 0;

}

#endif /* tests.hpp */