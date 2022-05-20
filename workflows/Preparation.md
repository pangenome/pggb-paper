```shell
mkdir -p ~/tools $$ cd ~/tools

#git clone --recursive https://github.com/waveygang/wfmash
cd wfmash
git checkout master && git pull && git submodule update --init --recursive
git checkout 48493fbe2d7e65b7e371b8940590c88dac13e4ad
cmake -H. -Bbuild && cmake --build build -- -j 48
mv build/bin/wfmash build/bin/wfmash-48493fbe2d7e65b7e371b8940590c88dac13e4ad
cd ..

#git clone --recursive https://github.com/ekg/seqwish.git
cd seqwish
git checkout master && git pull && git submodule update --init --recursive
git checkout f209482924bbb7763d927029ee63cd16c5ddb44c
cmake -H. -Bbuild && cmake --build build -- -j 48
mv bin/seqwish bin/seqwish-f209482924bbb7763d927029ee63cd16c5ddb44c
cd ..

#git clone --recursive https://github.com/pangenome/smoothxg.git
cd smoothxg
git checkout master && git pull && git submodule update --init --recursive
git checkout a2a15538f9e5b03ba2eb5a8f2d8d91f2ade28cd2
cmake -H. -Bbuild && cmake --build build -- -j 48
mv bin/smoothxg bin/smoothxg-a2a15538f9e5b03ba2eb5a8f2d8d91f2ade28cd2
cd ..

#git clone --recursive https://github.com/pangenome/odgi.git
cd odgi
git checkout master && git pull && git submodule update --init --recursive
git checkout 4d4acae6d7e27fed7cdbfc5c23c3e32d28092e89
cmake -H. -Bbuild && cmake --build build -- -j 48
mv bin/odgi bin/odgi-4d4acae6d7e27fed7cdbfc5c23c3e32d28092e89
cd ..

#git clone --recursive https://github.com/pangenome/pggb.git
cd pggb
git checkout master && git pull && git submodule update --init --recursive
git checkout 4f71d69f4c0bcacbd9b286c1aaed86e7c37fc0a2
sed 's,"$fmt" wfmash,"$fmt" ~/tools/wfmash/build/bin/wfmash-48493fbe2d7e65b7e371b8940590c88dac13e4ad,g' pggb -i
sed 's,"$fmt" seqwish,"$fmt" ~/tools/seqwish/bin/seqwish-f209482924bbb7763d927029ee63cd16c5ddb44c,g' pggb -i
sed 's,"$fmt" smoothxg,"$fmt" ~/tools/smoothxg/bin/smoothxg-a2a15538f9e5b03ba2eb5a8f2d8d91f2ade28cd2,g' pggb -i
sed 's,"$fmt" odgi,"$fmt" ~/tools/odgi/bin/odgi-4d4acae6d7e27fed7cdbfc5c23c3e32d28092e89,g' pggb -i
mv pggb pggb-4f71d69f4c0bcacbd9b286c1aaed86e7c37fc0a2
cd ..
```