#!/bin/env zsh -e

cd ..

echo 'SRAの名前を入力してください(.sraの前)'
read accession

echo '@edu.k.u-tokyo.ac.jpのアドレスを設定すると、プログラム終了時にメールが届きます。メールを設定しますか？[y/n]'
read mail
if [ $mail = 'y' ]; then
  echo '@edu.k.u-tokyo.ac.jpの前を入力してください'
  read address
else
  echo 'メールを設定しません'
fi

cd sra

prefetch ${accession} --max-size 100GB

mv ./${accession}/${accession}.sra .
rm -r ${accession}

echo 'ダウンロードが完了しました。data.csvに詳細を記入してください。'

#終了時メールで送信
if [ $mail = 'y' ]; then
  echo "Process done　メールを送信しました。" | mail -s "Process done" ${address}@edu.k.u-tokyo.ac.jp
else
  echo "Process done"
fi