# Developer Documentation

## Local Compilation

Build from source with Cargo with:

```bash
cargo build --release
./target/release/krak --help
```

## Local Testing

To ensure all tests pass, run:

```bash
cargo test --all --no-fail-fast --verbose
```

## Local Linting and Formatting

To check the format and lint of the code, run:

```bash
cargo fmt -- --check
cargo clippy --all-targets --all-features -- -D warnings
```
