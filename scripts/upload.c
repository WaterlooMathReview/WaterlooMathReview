#include <stdlib.h>
#include <unistd.h>

/*
int main() {
  execl("/users/mathreview/www/scripts/upload.py", "upload.py");
}*/

#define _POSIX_C_SOURCE 200112L // Needed with glibc (e.g., linux).
#include <stdio.h>
#include <sys/types.h>
#include <unistd.h>

void report (uid_t real) {
    printf (
        "Real UID: %d Effective UID: %d\n",
        real,
        geteuid()
    );
}

int main (void) {
    uid_t real = getuid();
    report(real);
    seteuid(real);
    report(real);
    return 0;
}
