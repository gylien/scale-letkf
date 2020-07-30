LIST_DIR = "dafcst_ncl"
NC_DIR = "dafcst_nc"
POSTFIX = ".nc"
SCRP_DIR = "/work/hp150019/share/saitama_aip_uploader"  
###DEST = "pawr-netcdf-saitama-now"  ### OPE
DEST = "pawr-netcdf-saitama-now-syd" ### TEST

operational = true

old = []
while true
  t = Time.now
  latest = Dir.glob("./#{LIST_DIR}/*").sort
  diff = latest - old
  if diff.length > 0
    diff.each{|f|
      next if /\.done$/ =~ f
      ff = File.join(NC_DIR, File.basename(f)) << POSTFIX
      puts "#{t}: #{ff}, #{File.mtime(ff)}"

      # TRANSFER DATA
#      system("scp -p -P 2222 #{ff} nowcast_pawr@hibiki.r-ccs.riken.jp:/global_hibiki/nowcast/") || raise 
      system("echo python #{SCRP_DIR}/aip.py --file #{ff} --bucket #{DEST}}") || raise 
      
      if operational
        File.unlink(f) # FOR OPERATIONAL MODE
      else
        File.rename(f, "#{f}.done") # FOR TEST MODE
      end
    }
  end
  old = latest
  sleep 1
end

