データ転送テスト用スクリプト

# まず、OFP にログインする前に、hibiki.r-ccs.riken.jp の鍵認証情報を ssh-add する。
# この部分は環境によるので、必要に応じて行う。
$ ssh-agent >> SSH_AGENT_SCRIPT
$ . SSH_AGENT_SCRIPT

# 事前に鍵をもっておく必要がある。
$ ssh-add ~/.ssh/id_rsa_nowcast_pawr_hibiki.aics.riken.jp

# 鍵情報を転送するモードでSSHする (-A オプション)。
$ ssh -A d41002@ofp03.jcahpc.jp

# 試験ディレクトリに移動。
$ cd /work/hp150019/share/otsuka/real-time-output-transfer

# 起動。
$ ruby transfer-fcst-test.rb

# データは hibiki.r-ccs.riken.jp:/global_hibiki/nowcast/ に送られる。
