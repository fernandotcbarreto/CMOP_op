FROM alpine:latest AS build

RUN apk add --no-cache \
    alpine-sdk \ 
    cmake \
    gfortran

COPY code /code

WORKDIR /build

RUN cmake /code && cmake --build . && cmake --install .

FROM alpine:latest

COPY --from=build /usr/local/bin/cmop /usr/bin/
COPY --from=build /usr/lib/libgcc_s.so.1 /usr/lib/
COPY --from=build /usr/lib/libgfortran.so.5 /usr/lib/
COPY --from=build /usr/lib/libquadmath.so.0 /usr/lib/
COPY --from=build /lib/ld-musl-x86_64.so.1 /lib/

WORKDIR /run

CMD ["cmop"]