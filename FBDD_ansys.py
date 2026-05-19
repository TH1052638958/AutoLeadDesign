prompt = 'Generate a novel valid molecule SMILES which contains one fragment of [ N[C@@H]1Cc2ccc3ncccc3c2[C@H](OC2C=Cc3ccccc32)C1 , Oc1ccc2c(c1)CNC2c1nc(N)c2ccccc2n1,C1C2=C3NC(C)=CC=C3C1CN2 ] at least and do not generate any English text.'

history_message=[{'role': 'user', 'content': prompt},
                 {'role': 'assistant', 'content': 'O=C(NC[C@@H]1Cc2ccc3ncccc3c2[C@H](OC2C=Cc3ccccc32)C1)N[C@H]1CCC[C@H]1c1ccc2ncccc2c1'},
                 {'role': 'user', 'content': 'Explain  the principles of fragment-based drug design that guided your decision-making in above process, including which fragments were chosen and how you design the compounds based on the fragments you choose.'}]


import os
from openai import OpenAI
import json
from openai.types.chat import completion_create_params

os.environ["OPENAI_BASE_URL"] = 'https://api.deepseek.com/v1'
os.environ["OPENAI_API_KEY"] = 'sk-6336b51c54f84a51b9262ba593f4997b'
model_engine = 'deepseek-chat'  # You can choose a different model if desired
client = OpenAI()
try:
    completion = client.chat.completions.create(
        model=model_engine,
        messages=history_message,
        n=1,
        # max_tokens=60,
        temperature=0.6,
        stop="!",
        user="user"
    )
except Exception as e:
    print(f'API ERROR:{e},retrying.......')

result = json.loads(completion.model_dump_json())["choices"][0]["message"]["content"]
print(result)