#include "PDF2DB.h"

void PDF2DB::Clear(Option_t *opt) {
  for (auto it : Pdf2Map)
    delete it.second;
  Pdf2Map.clear();
}

void PDF2DB::Print(Option_t *opt) const {
  for (auto it : Pdf2Map)
    printf("PDF2DB::Print %s\n", it.first.c_str());
}

void PDF2DB::Draw(Option_t *opt) {
  for (auto it : Pdf2Map)
    it.second->Draw();
}

void PDF2DB::Add(std::string name, PDF2 *pdf) { Pdf2Map[name] = pdf; }

PDF2 *PDF2DB::Get(std::string name) {
  std::map<std::string, PDF2 *>::iterator it = Pdf2Map.find(name);
  if (it == Pdf2Map.end())
    return 0;
  return it->second;
}
