#!/bin/env zsh -e

cd ..

echo 'data.csvに記入したSRAをダウンロードします。？[y/n]'
read qestion
if [ $qestion = 'y' ]; then
  echo 'プログラムを実行します'
else
  echo 'プログラムを終了します'
  exit
fi


echo '@edu.k.u-tokyo.ac.jpのアドレスを設定すると、プログラム終了時にメールが届きます。メールを設定しますか？[y/n]'
read mail
if [ $mail = 'y' ]; then
  echo '@edu.k.u-tokyo.ac.jpの前を入力してください'
  read address
else
  echo 'メールを設定しません'
fi

cut -d ',' -f 1 data.csv | sed -e '1d' | while read line

do
 find ./sra/${line}.sra
    if [ $? = 0 ] ; then
      echo 'SRAがあります。'
    else
      echo 'SRAをダウンロードします。'
      cd sra
      prefetch ${line} --max-size 100GB
      mv ./${line}/${line}.sra .
      rm -r ${line}
      cd ..
    fi
 
done

#終了時メールで送信
if [ $mail = 'y' ]; then
  echo "Process done　メールを送信しました。" | mail -s "Process done" ${address}@edu.k.u-tokyo.ac.jp
else
  echo "Process done"
fi

