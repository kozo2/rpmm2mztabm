# rpmm2mztabm

RIKEN Plant Metabolome Metadatabase のデータと解析成果物を管理します。

## ファイル構成

```text
projects/
  RPMM0003/
    data/                 # 取得した元データ（元の内部ディレクトリ構造を保持）
    tables/               # 定量表 Parquet と試料メタデータ TOML
    phenodata/            # このプロジェクトの列・値・セル注釈、要確認一覧、出典
    mztabm/               # mzTab-M を生成した場合の保存先
  RPMM0004/
  ...
reports/
  phenodata/              # 全プロジェクトの注釈・品質検証・集計
  tables/                 # 定量表の全件レポートと要確認一覧
  relocation/             # 旧→新配置、元のSHA256、再生成の照合結果
investigation/            # 再生成スクリプト、取得記録、辞書、調査資料
```

入力データは47プロジェクトにあり、定量表は43プロジェクト、今回対象の
`phenodata.csv`（大文字・小文字を区別しない）の注釈は31プロジェクトにあります。
各プロジェクトの成果物は、その入力データと同じフォルダに配置しています。

`RPMM0003/AnalysedDataset/phenodata.csv` のような記録内のパスは、配置に
依存しない論理キーです。実ファイルは
`projects/RPMM0003/data/AnalysedDataset/phenodata.csv` にあります。
論理キーを維持することで、セルID・対応表ID・出典参照を保っています。

## Phenodata 注釈の再生成

Python 3.9以降の標準ライブラリだけで実行でき、通常の生成では通信しません。
以下はリポジトリルートからのコマンドです。スクリプトのパスを絶対パスにすれば、
作業ディレクトリを変えても同じ入力を読み込みます。

```sh
python3 investigation/annotate_phenodata.py
PYTHONPATH=investigation python3 -m phenodata_ontology.validate_export
python3 -m unittest investigation.phenodata_ontology.test_annotation
python3 investigation/verify_relocation.py
```

全件集計は `reports/phenodata/`、プロジェクト別の同じ内容は
`projects/RPMMxxxx/phenodata/` に出力されます。分類ルール・オントロジーの
版・出典については [注釈処理の説明](investigation/phenodata_ontology/README.md) を参照してください。

## 定量表の再現性確認

[固定した依存バージョン](investigation/requirements-export.txt) をPython 3.12の
仮想環境へインストールして実行します。元成果物の日付は `2026-10-27` です。
以下は別出力先に再生成して既存成果物と比較する例です。

```sh
python3.12 -m venv /private/tmp/rpmm-repro-venv
/private/tmp/rpmm-repro-venv/bin/python -m pip install -r investigation/requirements-export.txt
/private/tmp/rpmm-repro-venv/bin/python investigation/export_tables.py \
  --output-root /private/tmp/rpmm-table-rebuild \
  --reports-root /private/tmp/rpmm-table-rebuild-reports \
  --generated-on 2026-10-27
python3 investigation/verify_relocation.py --rebuilt-tables-root /private/tmp/rpmm-table-rebuild
```

再配置時は106組のTOML・Parquet、計212ファイルを全件再生成し、移動前の
成果物とバイト単位で一致しました。元データのSHA256と、注釈全26,997セルの
原値・セルID・対応表IDの保持も検証しています。

[再配置の対応表](reports/relocation/manifest.json)、
[検証結果](reports/relocation/validation.json)、
[テーブルの再生成結果](reports/relocation/table_reproduction.json) を保存しています。
