import os
import torch
from torch import nn
from torch.utils.data import DataLoader, Dataset
from torch.nn import Softmax
import pickle as pkl
import argparse


class EsmClassifier(nn.Module):
    def __init__(
        self,
        embed_size=128,
    ):
        super(EsmClassifier, self).__init__()

        self.fc1 = nn.Linear(embed_size, embed_size)
        self.fc2 = nn.Linear(embed_size, 2)

    def forward(self, src):
        x = self.fc1(src)
        out = self.fc2(x)
        return out


class VirDataset(Dataset):
    def __init__(self, embed_file):
        self.inputs = torch.load(embed_file, weights_only=False)
        self.embed_size = self.inputs[0][1].shape[0]


    def __len__(self):
        return len(self.inputs)

    def __getitem__(self, index):
        return self.inputs[index][1], self.inputs[index][0]


parser = argparse.ArgumentParser()
parser.add_argument('--input', type=str)
parser.add_argument('--db', type=str)
parser.add_argument('--output', type=str)
args = parser.parse_args()
input = args.input
db = args.db
output = args.output

test_dataset = VirDataset(input)
embed_size = test_dataset.embed_size

test_loader = DataLoader(
    dataset=test_dataset,
    batch_size=16,
    shuffle=True,)

device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

model = EsmClassifier(
    embed_size=embed_size,
).to(device)

if torch.cuda.device_count() > 1:
    model = nn.DataParallel(model)
model.to(device)

checkpoint = torch.load(os.path.join(db, 'classifier.pt'))
model.load_state_dict(checkpoint['state_dict'])

softmax = Softmax(dim=0)
_ = model.eval()
result = {}

with torch.no_grad():
    for step, (batch_x, labels) in enumerate(test_loader):
        scores = model(batch_x.to(device)).cpu().numpy()

        for i in torch.arange(0, len(labels)):
            value = softmax(torch.tensor([scores[i][0], scores[i][1]])).tolist()
            label = labels[i]
            
            result[label] = value[1]
pkl.dump(result, open(os.path.join(output, 'protein_result.pkl'), 'wb'))
