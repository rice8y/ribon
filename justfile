wasm := "target/wasm32-unknown-unknown/release/ribon_plugin.wasm"
qa-dir := "target/qa"
typst-test-dir := qa-dir / "typst"
typage-docs-dir := "ribon-docs"

default:
  @just --list

plugin:
  cargo build --release --locked --target wasm32-unknown-unknown -p ribon-plugin
  cp {{wasm}} package/ribon_plugin.wasm

wasm-sync-test:
  ./scripts/check-wasm-sync.sh

wasm-test: plugin
  python3 scripts/validate-wasm.py

license-test:
  python3 scripts/validate-licenses.py

contrast-test:
  python3 scripts/validate-contrast.py

fmt-test:
  cargo fmt --all -- --check

test:
  cargo test --workspace --locked
  cargo clippy --workspace --all-targets --locked -- -D warnings

typst-test: plugin contrast-test
  mkdir -p {{typst-test-dir}}
  typst compile --root . tests/typst/smoke.typ {{typst-test-dir}}/smoke.pdf
  typst compile --root . tests/typst/api.typ {{typst-test-dir}}/api.pdf
  typst compile --root . tests/typst/api_all.typ {{typst-test-dir}}/api_all.pdf
  typst compile --root . tests/typst/layouts.typ {{typst-test-dir}}/layouts.pdf
  typst compile --root . tests/typst/extended.typ {{typst-test-dir}}/extended.pdf
  typst compile --root . tests/typst/publication.typ {{typst-test-dir}}/publication.pdf
  typst compile --root . tests/typst/plot_quality.typ {{typst-test-dir}}/plot-quality.pdf
  typst compile --root . tests/typst/plot_layout_quality.typ {{typst-test-dir}}/plot-layout-quality.pdf
  python3 scripts/validate-extended-pdf.py
  python3 scripts/validate-publication-render.py
  python3 scripts/validate-plot-quality.py
  python3 scripts/validate-plot-layout.py

plot-test: plugin
  mkdir -p {{typst-test-dir}}
  typst compile --root . tests/typst/plot_quality.typ {{typst-test-dir}}/plot-quality.pdf
  typst compile --root . tests/typst/plot_layout_quality.typ {{typst-test-dir}}/plot-layout-quality.pdf
  python3 scripts/validate-plot-quality.py
  python3 scripts/validate-plot-layout.py

plot-golden-update: plugin
  mkdir -p {{typst-test-dir}}
  typst compile --root . tests/typst/plot_quality.typ {{typst-test-dir}}/plot-quality.pdf
  typst compile --root . tests/typst/plot_layout_quality.typ {{typst-test-dir}}/plot-layout-quality.pdf
  python3 scripts/validate-plot-quality.py --update-golden
  python3 scripts/validate-plot-layout.py --update-golden

example: plugin
  mkdir -p {{qa-dir}}/examples
  typst compile --root . package/examples/secondary-structure.typ {{qa-dir}}/examples/secondary-structure.pdf

real-data-pdf: plugin
  mkdir -p {{qa-dir}}
  typst compile --root . tests/typst/real_data_render.typ {{qa-dir}}/ribon-real-data-validation.pdf

publication-pdf: plugin
  mkdir -p {{qa-dir}}
  typst compile --root . tests/typst/publication.typ {{qa-dir}}/ribon-publication-validation.pdf

render-test: real-data-pdf
  python3 scripts/validate-rendered-pdf.py
  python3 scripts/validate-render-golden.py

render-golden-update: real-data-pdf
  python3 scripts/validate-render-golden.py --update-golden

performance-test: plugin
  python3 scripts/validate-performance.py

pseudoknot-test: plugin
  python3 scripts/validate-pseudoknot.py

conditional-density2-test: plugin
  python3 scripts/validate-conditional-density2.py

conditional-density2-performance-test:
  python3 scripts/validate-conditional-density2-performance.py

exact-feature-test: plugin
  python3 scripts/validate-exact-features-real.py

release-check: fmt-test test wasm-sync-test wasm-test license-test typst-test pseudoknot-test conditional-density2-test conditional-density2-performance-test exact-feature-test performance-test render-test

docs:
  just --justfile package/justfile docs

docs-current:
  just --justfile package/justfile docs-current

images:
  just --justfile package/justfile images

install:
  just --justfile package/justfile install

typage-init:
  typage init {{typage-docs-dir}} \
    --starter github:rice8y/typage-starter-typst-docs \
    --install-plugins
  rm -rf "{{typage-docs-dir}}/package"
  mkdir -p "{{typage-docs-dir}}/package"
  cp package/LICENSE package/NOTICE.md package/README.md package/THIRD_PARTY.md package/justfile package/lib.typ package/ribon_plugin.wasm package/typst.toml "{{typage-docs-dir}}/package/"
  cp -R package/src "{{typage-docs-dir}}/package/"

typage-dev:
  typage dev --root {{typage-docs-dir}}

clean-qa:
  rm -rf {{qa-dir}}
